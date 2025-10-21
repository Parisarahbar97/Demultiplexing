# Demultiplexing & Genotype Imputation Workflow

This repository now captures the full genotype → imputation → demultiplexing
pipeline that we iterated on manually.  The workflow is split into three
logical phases:

1. **Pre-imputation QC (PLINK2 + bcftools)** – produces clean per-pool VCFs.
2. **MIS preparation / post-processing** – combines pools, prepares per-chr
   upload files, and merges MIS results after imputation.
3. **Demultiplexing** – existing Nextflow pipeline (`main.nf`) that runs
   `harmonize_vcf_bam`, `dsc_pileup`, `demuxlet`, and QC.

All heavy lifting happens inside the `Docker/Dockerfile.impute` image which
contains PLINK 2, PLINK 1.9, bcftools, samtools, unzip, pigz, gawk, etc.

## 0. Build the helper image

```bash
cd Docker
docker build -t parisa/genotype:impute -f Dockerfile.impute .
```

## 1. Pre-imputation QC

Input manifest: one VCF per line (absolute paths), e.g.
`imputation_work/00_raw/vcf_manifest.txt`.

```bash
docker run --rm -v /home/pr422:/host parisa/genotype:impute bash -lc '
  cd /host/RDS/live/Users/Parisa/Demultiplexing &&
  ./scripts/run_plink_qc.sh \
    --manifest /host/RDS/live/Users/Parisa/imputation_work/00_raw/vcf_manifest.txt \
    --outdir   /host/RDS/live/Users/Parisa/imputation_work/01_plink_qc \
    --logdir   /host/RDS/live/Users/Parisa/imputation_work/logs \
    --threads  8'
```

Outputs (`imputation_work/01_plink_qc/`):

- `POOL.pgen/.pvar/.psam` – PLINK 2 binary dataset.
- `POOL.qcmetrics.smiss/acount` – basic QC tables.
- `POOL.clean.vcf.gz` – autosomal, biallelic, MAC ≥ 1 SNPs (used downstream).

## 2. Prepare MIS upload set (GRCh38)

Combine all pools, keep a single copy of each donor, merge, and slice by chr.

```bash
docker run --rm -v /home/pr422:/host parisa/genotype:impute bash -lc '
  cd /host/RDS/live/Users/Parisa/Demultiplexing &&
  ./scripts/mis_prepare.sh \
    --input-dir  /host/RDS/live/Users/Parisa/imputation_work/01_plink_qc \
    --output-dir /host/RDS/live/Users/Parisa/imputation_work/02_mis_prep \
    --threads    4'
```

Outputs (`imputation_work/02_mis_prep/`):

- `unique_vcfs.list`, `vcf_unique_manifest.tsv` – provenance of retained donors.
- `all_donors.samples.txt` – combined donor list (64 in our case).
- `all_donors.hg38.chr.vcf.gz` – merged multi-sample VCF.
- `by_chrom/all_donors.chr{1..22}.vcf.gz` (+ `.tbi`) – upload these to the
  Michigan Imputation Server (GRCh38 panel, e.g. 1000g-phase3-low).

**Manual step:** upload the 22 `.vcf.gz` files (and their `.tbi`) via your MIS
account. Choose an hg38 panel (e.g. `1000g-phase3-low (hg38)`), Eagle v2.4
prephasing, and Minimac4.

## 3. Post-imputation QC (after MIS download)

Point the script at the folder containing `chr_*.zip` from MIS. The helper
converts each chromosome archive to BCF (for robust concatenation), merges
them, applies the INFO/R2 filter, and optionally enforces an AF window.

```bash
docker run --rm -v /home/pr422:/host parisa/genotype:impute bash -lc '
  cd /host/RDS/live/Users/Parisa/Demultiplexing &&
  ./scripts/mis_postprocess.sh \
    --input-dir  /host/RDS/live/Users/Parisa/imputation_work/03_imputed/mis_job1_raw \
    --output-dir /host/RDS/live/Users/Parisa/imputation_work/03_imputed/mis_job1_post \
    --prefix     job1 \
    --r2-min     0.4 \
    --maf-min    0.0 \
    --maf-max    1.0'
```

Outputs (example: `imputation_work/03_imputed/mis_job1_post/post/`):

- `job1.merged.bcf` – concatenated chr1..22 BCF.
- `job1.R2filt.vcf.gz` – INFO/R2 filtered autosomal SNPs (no AF window).
- `job1.R2filt.MAF.vcf.gz` – optional AF-filtered VCF (useful if you impose an
  AF window via `--maf-min/--maf-max`).

### Optional: subset per pool

Using the manifest created in step 2, generate four-donor VCFs:

```bash
POST=/home/pr422/RDS/live/Users/Parisa/imputation_work/03_imputed/mis_job1_post/post/job1.R2filt.vcf.gz
MAN=/home/pr422/RDS/live/Users/Parisa/imputation_work/02_mis_prep/vcf_unique_manifest.tsv
POOL_OUT=/home/pr422/RDS/live/Users/Parisa/imputation_work/03_imputed/mis_job1_pools
mkdir -p "${POOL_OUT}"
for pool in D10A D10P D11A D11P D12A D12P D13A D13P D1A D1P D2A D2P \
            D3A D3P D4A D4P D5A D5P D6A D6P D7A D7P D8A D8P D9A D9P; do
  awk -v pfx="${pool}.clean.vcf.gz" -F"\t" '$1==pfx{print $2}' "${MAN}" > "${POOL_OUT}/${pool}.donors.txt"
  if [[ -s "${POOL_OUT}/${pool}.donors.txt" ]]; then
    bcftools view -S "${POOL_OUT}/${pool}.donors.txt" "${POST}" -v snps -m2 -M2 -c 1:minor -Oz -o "${POOL_OUT}/${pool}.imputed.vcf.gz"
    tabix -f -p vcf "${POOL_OUT}/${pool}.imputed.vcf.gz"
  fi
done
```

Use these per-pool VCFs (GP preferred) in the Nextflow demultiplexing pipeline.

## 4. Demultiplexing (existing workflow)

Update `examples/samples.csv` to point to the imputed VCFs (GP/DS) and run:

```bash
nextflow run main.nf -profile dsi \
  --samples /home/pr422/RDS/live/Users/Parisa/Demultiplexing/examples/samples.csv \
  --outdir  /home/pr422/RDS/live/Users/Parisa/demux_out_nf \
  -with-report -with-timeline -with-trace
```

The pipeline stages are unchanged (`harmonize_vcf_bam`, `dsc_pileup`,
`run_demuxlet`, `qc_demuxlet`), but you can now feed high-quality imputed VCFs
with GP/DS instead of hard calls.

## Repo structure highlights

```
Docker/
  Dockerfile.impute      # PLINK2 + bcftools + helper utilities
examples/
  samples.csv            # Sample sheet used by main.nf
modules/
  *.nf                   # Nextflow DSL2 modules (harmonize, pileup, demuxlet, QC)
scripts/
  run_plink_qc.sh        # Step 1 – per-pool PLINK2 QC
  mis_prepare.sh         # Step 2 – combine pools & emit chr slices
  mis_postprocess.sh     # Step 3 – merge/filter MIS outputs
main.nf                  # Demultiplexing workflow
nextflow.config          # Profiles and container definitions
run_command.sh           # Ready-to-run command snippets
```

Feel free to tailor thresholds (`--r2-min`, `--maf-min`) or add extra QC
steps (e.g. HWE post-imputation) depending on the downstream analysis.
