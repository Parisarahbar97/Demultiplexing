#!/usr/bin/env bash
# Example command snippets for the full genotype → imputation → demultiplexing workflow.
# Adjust paths to fit your environment. All commands assume this repository as CWD.

set -euo pipefail

# 0) Change to repo root
cd /home/pr422/RDS/live/Users/Parisa/Demultiplexing

# 1) Pre-imputation QC (inside Docker)
docker run --rm -v /home/pr422:/host parisa/genotype:impute bash -lc '\
  cd /host/RDS/live/Users/Parisa/Demultiplexing && \
  ./scripts/run_plink_qc.sh \
    --manifest /host/RDS/live/Users/Parisa/imputation_work/00_raw/vcf_manifest.txt \
    --outdir   /host/RDS/live/Users/Parisa/imputation_work/01_plink_qc \
    --logdir   /host/RDS/live/Users/Parisa/imputation_work/logs \
    --threads  8'

# 2) Prepare TOPMed upload set (chr1..chr22 VCFs)
docker run --rm -v /home/pr422:/host parisa/genotype:impute bash -lc '\
  cd /host/RDS/live/Users/Parisa/Demultiplexing && \
  ./scripts/mis_prepare.sh \
    --input-dir  /host/RDS/live/Users/Parisa/imputation_work/01_plink_qc \
    --output-dir /host/RDS/live/Users/Parisa/imputation_work/02_mis_prep \
    --threads    4'

# (Manual) Upload /host/RDS/live/Users/Parisa/imputation_work/02_mis_prep/by_chrom/all_donors.chr*.vcf.gz (+ .tbi)
#               to the Michigan Imputation Server (GRCh38 panel).

# 3) Post-imputation QC after MIS download
docker run --rm -v /home/pr422:/host parisa/genotype:impute bash -lc '\
  cd /host/RDS/live/Users/Parisa/Demultiplexing && \
  ./scripts/mis_postprocess.sh \
    --input-dir  /host/RDS/live/Users/Parisa/imputation_work/03_imputed/mis_job1_raw \
    --output-dir /host/RDS/live/Users/Parisa/imputation_work/03_imputed/mis_job1_post \
    --prefix     job1 \
    --r2-min     0.4 \
    --maf-min    0.0 \
    --maf-max    1.0'

# 4) Rebuild per-pool VCFs directly from the merged MIS output (required for demuxlet)
docker run --rm -v /home/pr422:/host parisa/genotype:impute bash -lc '\
  cd /host/RDS/live/Users/Parisa/Demultiplexing && \
  ./scripts/mis_make_pool_vcfs.sh \
    --merged-vcf /host/RDS/live/Users/Parisa/imputation_work/03_imputed/mis_job1_post/post/job1.R2filt.vcf.gz \
    --lists-dir  /host/RDS/live/Users/Parisa/vcf_per_samplepool/lists \
    --output-dir /host/RDS/live/Users/Parisa/imputation_work/03_imputed/mis_job1_pools_fix'

# 5) Run demultiplexing nextflow pipeline
nextflow run main.nf -profile dsi \
  --samples "/home/pr422/RDS/live/Users/Parisa/Demultiplexing/examples/samples.csv" \
  --outdir  "/home/pr422/RDS/live/Users/Parisa/demux_out_nf" \
  -with-report -with-timeline -with-trace
