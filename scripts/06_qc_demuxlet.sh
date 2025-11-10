#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 06_qc_demuxlet.sh --best demuxlet.best --barcodes FILE --pool ID [options]

Prints QC summaries (rows_vs_barcodes, status counts, singlet posteriors,
top doublet pairs) and optionally writes:
  * donor_QC.tsv (per-donor singlet/doublet counts)
  * assignments table compatible with the eQTL pipeline

Options:
  --best FILE         demuxlet.best
  --barcodes FILE     whitelist (same as pileup)
  --pool ID           Pool ID used to prefix barcodes in assignments output
  --outdir DIR        Directory for donor_QC.tsv (default: alongside best file)
  --assignments FILE  Path to write Barcode/Individual_Assignment/DropletType
  -h, --help          Show usage
USAGE
}

BEST=""
BARCODES=""
POOL=""
OUTDIR=""
ASSIGN=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --best) BEST="$2"; shift 2 ;;
    --barcodes) BAR CODES="$2"; shift 2 ;;
  esac
end while
