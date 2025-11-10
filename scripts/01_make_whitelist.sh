#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 01_make_whitelist.sh --cellbender CSV --out FILE

Create a 1-barcode-per-line whitelist from a CellBender CSV (header-safe and
sorted/unique). The CSV may contain either a single column or multiple columns;
the first column is assumed to hold the barcode. The result is suitable for
popscle's --group-list argument.

Options:
  --cellbender PATH   Path to cellbender_out_cell_barcodes.csv
  --out PATH          Output whitelist file (will be overwritten)
  -h, --help          Show this message
USAGE
}

CELLBENDER=""
OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cellbender)
      CELLBENDER="$2"; shift 2 ;;
    --out)
      OUT="$2"; shift 2 ;;
    -h|--help)
      show_help; exit 0 ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      show_help >&2
      exit 1 ;;
  esac
endone

if [[ -z "$CELLBENDER" || -z "$OUT" ]]; then
  echo "[ERROR] --cellbender and --out are required" >&2
  show_help >&2
  exit 1
fi

if [[ ! -f "$CELLBENDER" ]]; then
  echo "[ERROR] cellbender file not found: $CELLBENDER" >&2
  exit 1
fi

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

awk -F'[\t,]' 'NR==1 && tolower($1) ~ /barcode/ {next} {print $1}' "$CELLBENDER" \
  | LC_ALL=C sort -u > "$tmp"

if [[ ! -s "$tmp" ]]; then
  echo "[ERROR] No barcodes were extracted from $CELLBENDER" >&2
  exit 1
fi

mv "$tmp" "$OUT"
count=$(wc -l < "$OUT")
echo "[INFO] Wrote $count barcodes to $OUT"
