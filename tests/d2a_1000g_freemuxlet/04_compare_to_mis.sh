#!/usr/bin/env bash
# Compare the 1000G freemuxlet runs to the existing MIS-based freemux K=4 results.

set -euo pipefail

HOST=/home/pr422
OUTDIR=$HOST/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test
MIS_FREEMUX=$HOST/RDS/live/Users/Parisa/demux_manual/D2A/misVCF_run/freemux_K4.clust1.samples.gz
OUTTXT=$OUTDIR/compare_1000G_vs_MIS.txt

{
  echo "=== MIS-based freemux (K=4) ==="
  if [ -f "$MIS_FREEMUX" ]; then
    gzip -cd "$MIS_FREEMUX" | awk -F'\t' 'NR>1{c[$5]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' | sort
  else
    echo "(Missing: $MIS_FREEMUX)"
  fi
  echo
  for K in 3 4 5; do
    F="$OUTDIR/freemux1000G_K${K}.clust1.samples.gz"
    echo "=== 1000G-based freemux (K=$K) ==="
    if [ -f "$F" ]; then
      gzip -cd "$F" | awk -F'\t' 'NR>1{c[$5]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' | sort
    else
      echo "(Missing: $F)"
    fi
    echo
  done
} | tee "$OUTTXT"

echo "Comparison written: $OUTTXT"
