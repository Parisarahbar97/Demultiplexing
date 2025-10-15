#!/usr/bin/env bash
set -euo pipefail

OUTDIR=/home/pr422/RDS/live/Users/Parisa/demux_out

nextflow run main.nf -profile docker \
  --samples examples/samples.csv \
  --outdir "$OUTDIR" \
  --field GT \
  --doublet_prior 0.10 \
  --vcf_filter biallelic_poly \
  --pileup_vcf clean \
  --demux_vcf clean
