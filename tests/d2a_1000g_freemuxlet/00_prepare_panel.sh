#!/usr/bin/env bash
# Prepare a bgzipped + indexed 1000G VCF copy for the D2A test run.

set -euo pipefail

HOST=/home/pr422
OUTDIR=$HOST/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test
PANEL=$HOST/RDS/live/Users/Parisa/alex_output/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf
OUTVCF=$OUTDIR/1000G_sites.vcf.gz

mkdir -p "$OUTDIR"

docker run --rm -u "$(id -u)":"$(id -g)" -v "$HOST":"$HOST" parisa/genotype:2.6 bash -lc '
  set -euo pipefail
  OUTVCF="'$OUTVCF'"; PANEL="'$PANEL'"
  if [ ! -f "$PANEL" ]; then
    echo "Missing 1000G VCF: $PANEL" >&2
    exit 1
  fi
  if ! bcftools view -h "$PANEL" | grep -q "^##contig=<ID=chr"; then
    echo "WARN: 1000G VCF header does not show chr-style contigs; BAM uses chr*." >&2
    echo "      Adjust contig names before running pileup if needed." >&2
  fi
  bgzip -c "$PANEL" > "$OUTVCF"
  tabix -f -p vcf "$OUTVCF"
  bcftools index -f "$OUTVCF"
  echo "Prepared: $OUTVCF"
'
