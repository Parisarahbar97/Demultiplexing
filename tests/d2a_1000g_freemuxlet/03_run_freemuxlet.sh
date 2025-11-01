#!/usr/bin/env bash
# Run freemuxlet on the 1000G-based pileup for a small grid of K values.

set -euo pipefail

HOST=/home/pr422
OUTDIR=$HOST/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test
PLP=$OUTDIR/D2A_1000G_pileup
BARCODES=$OUTDIR/D2A.barcodes.txt

for f in "${PLP}.plp.gz" "${PLP}.var.gz" "${PLP}.umi.gz" "${PLP}.cel.gz" "$BARCODES"; do
  [ -f "$f" ] || { echo "Missing pileup component or barcode list: $f" >&2; exit 1; }
done

for K in 3 4 5; do
  docker run --rm -u "$(id -u)":"$(id -g)" -v "$HOST":"$HOST" \
    --entrypoint /opt/conda/bin/popscle parisa/demux:2.1 freemuxlet \
      --plp "$PLP" \
      --nsample "$K" \
      --group-list "$BARCODES" \
      --out "$OUTDIR/freemux1000G_K${K}"
done

{
  for K in 3 4 5; do
    F="$OUTDIR/freemux1000G_K${K}.clust1.samples.gz"
    if [ -f "$F" ]; then
      echo "==== 1000G freemux K=$K ===="
      gzip -cd "$F" | awk -F'\t' 'NR>1{c[$5]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' | sort
    fi
  done
} | tee "$OUTDIR/1000G_freemux_summaries.txt"

echo "Summaries written: $OUTDIR/1000G_freemux_summaries.txt"
