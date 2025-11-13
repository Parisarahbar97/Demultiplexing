#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 04_run_freemuxlet.sh --plp PREFIX --barcodes FILE --outdir DIR [options]

Runs freemuxlet for a list of K values and prints singlet/doublet counts.

Options:
  --plp PREFIX       Prefix from dsc-pileup (without .plp.gz suffix)
  --barcodes FILE    Whitelist
  --outdir DIR       Output directory
  --ks LIST          Comma-separated K values (default: 3,4,5)
  --host-root PATH   Host path for Docker bind (default: /home/pr422)
  --image NAME       popscle image (default: parisa/demux:2.1)
  -h, --help         Show usage
USAGE
}

PLP=""
BARCODES=""
OUTDIR=""
KS="3,4,5"
HOST_ROOT=${HOST_ROOT:-/home/pr422}
IMAGE=${POPSCLE_IMAGE:-parisa/demux:2.1}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plp) PLP="$2"; shift 2 ;;
    --barcodes) BARCODES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --ks) KS="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

for var in PLP BARCODES OUTDIR; do
  if [[ -z "${!var}" ]]; then
    echo "[ERROR] --${var,,} is required" >&2
    show_help >&2
    exit 1
  fi
done

mkdir -p "$OUTDIR"
IFS=',' read -r -a KLIST <<< "$KS"

for K in "${KLIST[@]}"; do
  OUT="$OUTDIR/freemux_K${K}"
  set -x
  docker run --rm -u $(id -u):$(id -g) \
    -v "$HOST_ROOT":"$HOST_ROOT" \
    --entrypoint /opt/conda/bin/popscle "$IMAGE" freemuxlet \
      --plp "$PLP" \
      --nsample "$K" \
      --group-list "$BARCODES" \
      --out "$OUT"
  set +x
  echo "==== freemux K=$K ===="
  if [[ -f "${OUT}.clust1.samples.gz" ]]; then
    gzip -cd "${OUT}.clust1.samples.gz" | awk -F'\t' 'NR>1{c[$5]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' | sort
  else
    echo "(missing ${OUT}.clust1.samples.gz)"
  fi
  echo
  unset c
endone
