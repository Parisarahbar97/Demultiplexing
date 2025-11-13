#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 03_qc_demuxalot.sh --assignments FILE [options]

Summarize Demuxalot assignments and optionally emit a pool-prefixed
3-column metadata file (Barcode, Individual_Assignment, Demuxalot_droplet_type).

Options:
  --assignments FILE  assignments.tsv(.gz) or assignments_refined.tsv(.gz)
  --pool ID           Pool ID to prefix barcodes (optional)
  --out FILE          Path to write the 3-column TSV (optional)
  --top INT           Number of top doublet combos to print (default: 5)
  -h, --help          Show usage
USAGE
}

ASSIGN=""
POOL=""
OUT=""
TOP=5

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assignments) ASSIGN="$2"; shift 2 ;;
    --pool) POOL="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --top) TOP="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

if [[ -z "$ASSIGN" ]]; then
  echo "[ERROR] --assignments is required" >&2
  show_help >&2
  exit 1
fi

if [[ ! -f "$ASSIGN" ]]; then
  echo "[ERROR] Missing assignments file: $ASSIGN" >&2
  exit 1
fi

read_data() {
  if [[ "$ASSIGN" == *.gz ]]; then
    gzip -dc "$ASSIGN"
  else
    cat "$ASSIGN"
  fi
}

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
read_data > "$tmp"

TOTAL=$(( $(wc -l < "$tmp") - 1 ))
if (( TOTAL <= 0 )); then
  echo "[ERROR] $ASSIGN appears empty" >&2
  exit 1
fi

echo "rows $TOTAL"

awk -F'\t' -v top="$TOP" '
  NR>1{
    assign=$2
    droplet="SNG"
    if (assign=="doublet" || assign=="Doublet" || index(assign,"+")>0) {
      droplet="DBL"
      dbl_total++
      dbl_combo[assign]++
    } else {
      sng_total++
      sng_counts[assign]++
    }
  }
  END{
    printf "SNG\t%d\n", sng_total+0
    printf "DBL\t%d\n", dbl_total+0
    print ""
    print "singlet_counts"
    for (d in sng_counts) printf "%s\t%d\n", d, sng_counts[d]
    print ""
    print "top_doublet_counts (top=" top ")"
    if (length(dbl_combo)==0) next
    PROCINFO["sorted_in"]="@val_num_desc"
    count=0
    for (combo in dbl_combo) {
      printf "%s\t%d\n", combo, dbl_combo[combo]
      if (++count>=top) break
    }
  }
' "$tmp"

if [[ -n "$OUT" ]]; then
  mkdir -p "$(dirname "$OUT")"
  {
    printf "Barcode\tIndividual_Assignment\tDemuxalot_droplet_type\n"
    awk -F'\t' -v pool="$POOL" 'NR>1{
      assign=$2
      droplet="SNG"
      if (assign=="doublet" || assign=="Doublet" || index(assign,"+")>0) {
        droplet="DBL"
      }
      barcode=$1
      if (pool!="") barcode=pool"_"barcode
      print barcode"\t"assign"\t" (droplet=="SNG"?"SNG":"DBL")
    }' "$tmp"
  } > "$OUT"
  echo "[INFO] Wrote summarized assignments to $OUT"
fi
