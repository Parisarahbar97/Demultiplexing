#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 06_qc_demuxlet.sh --best FILE --barcodes FILE --pool ID [options]

Print QC summaries for demuxlet.best and optionally write donor_QC.tsv and an
assignments table (Barcode, Individual_Assignment, Demuxlet_droplet_type).

Options:
  --best FILE         Path to demuxlet.best
  --barcodes FILE     Whitelist used for pileup
  --pool ID           Pool ID (prefix added to barcodes in assignments)
  --outdir DIR        Directory for donor_QC.tsv (default: dirname(best))
  --assignments FILE  Path to write the assignments table (optional)
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
    --barcodes) BARCODES="$2"; shift 2 ;;
    --pool) POOL="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --assignments) ASSIGN="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

for var in BEST BARCODES POOL; do
  if [[ -z "${!var}" ]]; then
    echo "[ERROR] --${var,,} is required" >&2
    show_help >&2
    exit 1
  fi
done

if [[ ! -f "$BEST" || ! -f "$BARCODES" ]]; then
  echo "[ERROR] Missing best or barcode file" >&2
  exit 1
fi

: "${OUTDIR:=$(dirname "$BEST")}" # default

TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT
if head -n1 "$BARCODES" | grep -qi 'barcode'; then
  tail -n +2 "$BARCODES" > "$TMP"
else
  cat "$BARCODES" > "$TMP"
fi

BARC_COUNT=$(wc -l < "$TMP")
BEST_ROWS=$(($(wc -l < "$BEST") - 1))
OVERLAP=$(comm -12 <(tail -n +2 "$BEST" | cut -f2 | sort) <(sort "$TMP") | wc -l)

echo "rows_vs_barcodes"
echo "barcodes $BARC_COUNT"
echo "best_rows $BEST_ROWS"
echo "overlap   $OVERLAP"
echo

echo "status_counts"
tail -n +2 "$BEST" | awk -F'\t' '{c[$5]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' | sort

echo
echo "sng_posteriors"
tail -n +2 "$BEST" | awk -F'\t' '$5=="SNG"{t++; if($12>=0.99) p99++; else if($12>=0.9) p90++; else if($12>=0.5) p50++; else p0++;}
  END{
    if(t>0) printf "SNG=%d  >=0.99=%d (%.1f%%)  0.9-<0.99=%d (%.1f%%)  0.5-<0.9=%d (%.1f%%)  <0.5=%d (%.1f%%)\n", t,p99,100*p99/t,p90,100*p90/t,p50,100*p50/t,p0,100*p0/t;
    else print "SNG=0"
  }'
echo
echo "top_dbl_pairs"
tail -n +2 "$BEST" | awk -F'\t' '$5=="DBL"{split($6,a,","); pair=(a[1]<a[2]?a[1]"+"a[2]:a[2]"+"a[1]); c[pair]++}
  END{for(k in c) printf "%d %s\n",c[k],k}' | sort -k1,1nr | head

echo
echo "donor_QC"
tail -n +2 "$BEST" | awk -F'\t' '
  $5=="SNG"{split($6,a,","); s[a[1]]++; seen[a[1]]=1}
  $5=="DBL"{split($6,a,","); d1[a[1]]++; d2[a[2]]++; seen[a[1]]=1; seen[a[2]]=1}
END{
  printf "donor\tsinglets\tdoublets_as_d1\tdoublets_as_d2\tdoublets_total\n"
  for(d in seen){sd=s[d]+0; x=d1[d]+0; y=d2[d]+0; printf "%s\t%d\t%d\t%d\t%d\n", d, sd, x, y, x+y}
}' | sort -k1,1 > "$OUTDIR/donor_QC.tsv"
cat "$OUTDIR/donor_QC.tsv"

echo "[INFO] donor QC written to $OUTDIR/donor_QC.tsv"

if [[ -n "$ASSIGN" ]]; then
  mkdir -p "$(dirname "$ASSIGN")"
  {
    echo '"Barcode" "Individual_Assignment" "Demuxlet_droplet_type"'
    tail -n +2 "$BEST" | awk -F'\t' -v pool="$POOL" '{
      split($6,a,",");
      donor = ($5=="DBL") ? a[1] "/" a[2] : a[1];
      printf "\"%s_%s\" \"%s\" \"%s\"\n", pool, $2, donor, $5
    }'
  } > "$ASSIGN"
  echo "[INFO] assignments written to $ASSIGN"
fi
