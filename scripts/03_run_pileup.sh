#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 03_run_pileup.sh --sample ID --bam BAM --vcf VCF --barcodes FILE --outdir DIR [options]

Runs popscle dsc-pileup for one pool (defaults: CB/UB tags, MQ 20, BQ 13, cap 40).

Options:
  --sample ID         Pool identifier (used in messages)
  --bam PATH          Coordinate-sorted BAM with .bai index
  --vcf PATH          BAM-ordered VCF/BCF (bgzip + .tbi/.csi)
  --barcodes FILE     Whitelist produced by 01_make_whitelist.sh
  --outdir DIR        Output directory (created if missing)
  --threads INT       OpenMP threads (default: 8)
  --tag-group STR     SAM tag for barcodes (default: CB)
  --tag-umi STR       SAM tag for UMIs (default: UB)
  --min-mq INT        Minimum mapping quality (default: 20)
  --min-bq INT        Minimum base quality (default: 13)
  --cap-bq INT        Base quality cap (default: 40)
  --host-root PATH    Host root to mount inside Docker (default: /home/pr422)
  --image NAME        popscle Docker image (default: parisa/demux:2.1)
  -h, --help          Show usage
USAGE
}

SAMPLE=""
BAM=""
VCF=""
BARCODES=""
OUTDIR=""
THREADS=8
TAG_GROUP=CB
TAG_UMI=UB
MIN_MQ=20
MIN_BQ=13
CAP_BQ=40
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
    --tag-group) TAG_GROUP="$2"; shift 2 ;;
    --tag-umi) TAG_UMI="$2"; shift 2 ;;
    --min-mq) MIN_MQ="$2"; shift 2 ;;
    --min-bq) MIN_BQ="$2"; shift 2 ;;
    --cap-bq) CAP_BQ="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

for var in SAMPLE BAM VCF BARCODES OUTDIR; do
  if [[ -z "${!var}" ]]; then
    echo "[ERROR] --${var,,} is required" >&2
    show_help >&2
    exit 1
  fi
done

for f in "$BAM" "$BAM.bai" "$VCF" "$BARCODES"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing file: $f" >&2
    exit 1
  fi
done

if [[ -f "$VCF.tbi" ]]; then
  VCF_IDX="$VCF.tbi"
elif [[ -f "$VCF.csi" ]]; then
  VCF_IDX="$VCF.csi"
else
  echo "[ERROR] Missing VCF index (.tbi or .csi) for $VCF" >&2
  exit 1
fi

mkdir -p "$OUTDIR"
OUTPFX="$OUTDIR/${SAMPLE}_pileup"

set -x
# shellcheck disable=SC2086
docker run --rm -u "$(id -u)":"$(id -g)" \
  -v "$HOST_ROOT":"$HOST_ROOT" \
  -e OMP_NUM_THREADS="$THREADS" \
  --entrypoint /opt/conda/bin/popscle "$IMAGE" dsc-pileup \
    --sam        "$BAM" \
    --vcf        "$VCF" \
    --group-list "$BARCODES" \
    --tag-group  "$TAG_GROUP" \
    --tag-UMI    "$TAG_UMI" \
    --min-MQ     "$MIN_MQ" \
    --min-BQ     "$MIN_BQ" \
    --cap-BQ     "$CAP_BQ" \
    --out        "$OUTPFX"
set +x

echo "[INFO] pileup complete: ${OUTPFX}.{plp,var,umi,cel}.gz"
