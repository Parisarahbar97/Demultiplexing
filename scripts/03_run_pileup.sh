#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 03_run_pileup.sh --sample ID --bam BAM --vcf VCF --barcodes FILE --outdir DIR [options]

Runs popscle dsc-pileup for one pool (CB/UB tags, MQ 30, BQ 20, cap 40).

Options:
  --sample ID         Pool identifier (used in log messages)
  --bam PATH          Path to Cell Ranger BAM
  --vcf PATH          BAM-ordered VCF/BCF (with GP/PL and index)
  --barcodes FILE     Whitelist from 01_make_whitelist.sh
  --outdir DIR        Output directory (will be created)
  --threads INT       OpenMP threads (default: 40)
  --host-root PATH    Host path to mount inside Docker (default: /home/pr422)
  --image NAME        popscle Docker image (default: parisa/demux:2.1)
  -h, --help          Show usage
USAGE
}

SAMPLE=""
BAM=""
VCF=""
BARCODES=""
OUTDIR=""
THREADS=40
HOST_ROOT=${HOST_ROOT:-/home/pr422}
IMAGE=${POPSCLE_IMAGE:-parisa/demux:2.1}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --bam) BAM="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --barcodes) BARCODES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
endone

for var in SAMPLE BAM VCF BARCODES OUTDIR; do
  if [[ -z "${!var}" ]]; then
    echo "[ERROR] --${var,,} is required" >&2
    show_help >&2
    exit 1
  fi
done

for f in "$BAM" "$BAM".bai "$VCF" "$VCF".tbi "$BARCODES"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing file: $f" >&2
    exit 1
  fi
done

mkdir -p "$OUTDIR"
OUTPFX="$OUTDIR/${SAMPLE}_pileup"

set -x
# shellcheck disable=SC2086
docker run --rm -u $(id -u):$(id -g) \
  -v "$HOST_ROOT":"$HOST_ROOT" \
  -e OMP_NUM_THREADS="$THREADS" \
  --entrypoint /opt/conda/bin/popscle "$IMAGE" dsc-pileup \
    --sam        "$BAM" \
    --vcf        "$VCF" \
    --group-list "$BARCODES" \
    --tag-group  CB \
    --tag-UMI    UB \
    --min-MQ     30 \
    --min-BQ     20 \
    --cap-BQ     40 \
    --out        "$OUTPFX"
set +x

echo "[INFO] pileup complete: ${OUTPFX}.{plp,var,umi,cel}.gz"
