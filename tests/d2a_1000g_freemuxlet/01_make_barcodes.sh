#!/usr/bin/env bash
# Generate the D2A barcode whitelist from CellBender output.

set -euo pipefail

HOST=/home/pr422
OUTDIR=$HOST/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test
CSV=$HOST/RDS/live/Users/Parisa/EPILEP/healthy/qc/output_latest/D2A/D2A_cellbender_output/cellbender_out_cell_barcodes.csv
OUT=$OUTDIR/D2A.barcodes.txt

mkdir -p "$OUTDIR"

if [ ! -f "$CSV" ]; then
  echo "Missing CellBender CSV: $CSV" >&2
  exit 1
fi

awk -F'[\,\t]' 'NR==1 && tolower($1) ~ /barcode/ {next} {print $1}' "$CSV" \
  | LC_ALL=C sort -u > "$OUT"

echo "Wrote barcodes: $OUT ("$(wc -l < "$OUT")" lines)"
