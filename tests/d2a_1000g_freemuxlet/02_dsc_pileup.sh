#!/usr/bin/env bash
# Run dsc-pileup for D2A using the prepared 1000G site panel.

set -euo pipefail

HOST=/home/pr422
OUTDIR=$HOST/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test
BAM=$HOST/RDS/live/Users/Parisa/EPILEP/healthy/mapping/output/D2A_mapped/outs/possorted_genome_bam.bam
VCF=$OUTDIR/1000G_sites.vcf.gz
BARCODES=$OUTDIR/D2A.barcodes.txt
OUTPFX=$OUTDIR/D2A_1000G_pileup

for f in "$BAM" "$VCF" "$BARCODES"; do
  [ -f "$f" ] || { echo "Missing input: $f" >&2; exit 1; }
done

docker run --rm -u "$(id -u)":"$(id -g)" \
  -v "$HOST":"$HOST" \
  --entrypoint /opt/conda/bin/popscle \
  parisa/demux:2.1 dsc-pileup \
    --sam        "$BAM" \
    --vcf        "$VCF" \
    --group-list "$BARCODES" \
    --tag-group  CB \
    --tag-UMI    UB \
    --min-MQ     30 \
    --min-BQ     20 \
    --cap-BQ     40 \
    --out        "$OUTPFX"

echo "Pileup written: ${OUTPFX}.{plp,var,umi,cel}.gz"
