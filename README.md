# Manual Demultiplexing Toolkit

This repository now tracks the lean, manual workflow we settled on while
processing the epileptic snRNA-seq pools. The pipeline is intentionally simple
(no Nextflow) and mirrors the steps we validated on D2A/D13A:

1. build a trustworthy barcode whitelist (CellBender → 1 line per barcode)
2. ensure the donor VCF uses the same contig order as the BAM
3. run `dsc-pileup` (CB/UB, MQ ≥30, BQ ≥20, cap 40) to collect SNP counts
4. run freemuxlet as a sanity check (discover true donor count K, DBL rate)
5. run demuxlet with the correct donor panel to obtain named assignments
6. generate QC summaries + the final metadata file (Barcode ↔ Individual_ID)

All heavy lifting happens inside two Docker images that you already use:

- `parisa/demux:2.1` – popscle (dsc-pileup / freemuxlet / demuxlet)
- `parisa/genotype:2.6` – bcftools, samtools, awk helpers

Set `HOST_ROOT=/home/pr422` (or export it once) so the scripts know how to
mount the filesystem inside Docker.

## Scripts

| Script | Purpose |
| --- | --- |
| `scripts/01_make_whitelist.sh` | header-safe whitelist from `cellbender_out_cell_barcodes.csv` |
| `scripts/02_reheader_vcf.sh` | reorder a pool VCF to match the BAM header |
| `scripts/03_run_pileup.sh` | run popscle `dsc-pileup` with CB/UB tags |
| `scripts/04_run_freemuxlet.sh` | run freemuxlet for `K=3,4,5` (configurable) |
| `scripts/05_run_demuxlet.sh` | run popscle `demuxlet` with realistic priors |
| `scripts/06_qc_demuxlet.sh` | produce QC summaries + eQTL-ready assignments |

Each script has `--help` describing the arguments. They are designed to be run
independently so you can stop, inspect outputs, and iterate before moving on.

## Quick-start (S1A example)

```bash
POOL=S1A
HOST_ROOT=/home/pr422
WD=$HOST_ROOT/RDS/live/Users/Parisa/demux_manual/$POOL/run1
mkdir -p "$WD"

# 1. whitelist
a=$HOST_ROOT/RDS/live/Users/Parisa/EPILEP/diseased/qc/output_latest/$POOL/${POOL}_cellbender_output/cellbender_out_cell_barcodes.csv
scripts/01_make_whitelist.sh --cellbender "$a" --out "$WD/${POOL}.barcodes.txt"

# 2. reheader VCF to BAM order (writes ${POOL}.bamorder.vcf.gz)
v=$HOST_ROOT/RDS/live/Users/Parisa/alex_output/epilep_cellranger_outputs/demuxlet_IO/genotype_inputs/${POOL}.vcf.gz
b=$HOST_ROOT/RDS/live/Users/Parisa/alex_output/epilep_cellranger_outputs/${POOL}_mapped/outs/possorted_genome_bam.bam
scripts/02_reheader_vcf.sh --vcf "$v" --bam "$b" --out "$WD/${POOL}.bamorder.vcf.gz"

# 3. pileup
scripts/03_run_pileup.sh \
  --sample "$POOL" \
  --bam "$b" \
  --vcf "$WD/${POOL}.bamorder.vcf.gz" \
  --barcodes "$WD/${POOL}.barcodes.txt" \
  --outdir "$WD"

# 4. freemuxlet sanity check
scripts/04_run_freemuxlet.sh \
  --plp "$WD/${POOL}_pileup" \
  --barcodes "$WD/${POOL}.barcodes.txt" \
  --outdir "$WD"
# Inspect singlet/doublet counts; decide on the correct K and roster.

# 5. demuxlet (adjust priors/filters as needed)
scripts/05_run_demuxlet.sh \
  --plp "$WD/${POOL}_pileup" \
  --vcf "$WD/${POOL}.bamorder.vcf.gz" \
  --barcodes "$WD/${POOL}.barcodes.txt" \
  --outdir "$WD/${POOL}_demuxlet" \
  --doublet-prior 0.05 \
  --min-total 80 --min-umi 40 --min-snp 30

# 6. QC + final metadata (Barcode, Individual_Assignment, droplet type)
scripts/06_qc_demuxlet.sh \
  --best "$WD/${POOL}_demuxlet/demuxlet.best" \
  --barcodes "$WD/${POOL}.barcodes.txt" \
  --pool "$POOL" \
  --outdir "$WD/${POOL}_demuxlet" \
  --assignments "$WD/${POOL}_demuxlet/${POOL}_assignments.tsv"

# Append assignments to the master file when satisfied:
#   cat ${POOL}_assignments.tsv >> /home/pr422/RDS/live/Users/Parisa/parisa_eqtl/demultiplexed_assignments/demuxlet_assignments.txt
```

## Notes & best practices

- **Whitelist must be non-empty.** If `01_make_whitelist.sh` reports zero
  barcodes, stop and inspect the CellBender output.
- **Keep the full imputed SNP set.** Do not filter sites with AC>0 across the
  small pool—doing so previously halved useful evidence and made demuxlet label
  most cells as doublets. Only reheader; do not drop variants unless you later
  create a polymorphic-only VCF on purpose.
- **Run freemuxlet on every pool.** It’s the fastest way to discover missing
  donors, wrong K, or poor pileups. If freemuxlet and demuxlet disagree, fix
  the roster (subset the merged MIS VCF to the donors freemuxlet found) and
  rerun.
- **UMI tag:** popscle expects UB for 10x data. Use `--tag-UMI UR` only if
  logs show “Cannot find UMI tag UB”.
- **Assignments format:** `06_qc_demuxlet.sh --assignments` produces the exact
  column order used by the eQTL workflow (`"Barcode" "Individual_Assignment"
  "Demuxlet_droplet_type"`). Load this table into Seurat metadata to add the
  `Individual_ID` column.

This repo intentionally focuses on the demultiplexing steps we actually use
now. The older PLINK/MIS helpers and Nextflow workflow were removed so the
instructions stay small and repeatable.
