#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 05_run_demuxlet.sh --plp PREFIX --vcf VCF --barcodes FILE --outdir DIR [options]

Runs popscle demuxlet with configurable priors and cell-level filters.

Required:
  --plp PREFIX       Output prefix from dsc-pileup (same as used by freemuxlet)
  --vcf PATH         BAM-ordered VCF (bgzip + index)
  --barcodes FILE    Whitelist
  --outdir DIR       Output directory (demuxlet.best etc.)

Optional:
  --sm-list FILE     File listing donor IDs (default: derived from VCF)
  --doublet-prior F  Prior doublet rate (default: 0.05)
  --min-total INT    --min-total passed to demuxlet (default: 0)
  --min-umi INT      --min-umi (default: 0)
  --min-snp INT      --min-snp (default: 0)
  --field NAME       FORMAT field (GP, DS, PL; default: GP)
  --host-root PATH   Host path for Docker bind (default: /home/pr422)
  --image NAME       popscle image (default: parisa/demux:2.1)
  -h, --help         Show usage
USAGE
}

PLP=""
VCF=""
BARCODES=""
OUTDIR=""
SM_LIST=""
DOUBLETP=0.05
MIN_TOTAL=0
MIN_UMI=0
MIN_SNP=0
FIELD=GP
HOST_ROOT=${HOST_ROOT:-/home/pr422}
IMAGE=${POPSCLE_IMAGE:-parisa/demux:2.1}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plp) PLP="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --barcodes) BARCODES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --sm-list) SM_LIST="$2"; shift 2 ;;
    --doublet-prior) DOUBLETP="$2"; shift 2 ;;
    --min-total) MIN_TOTAL="$2"; shift 2 ;;
    --min-umi) MIN_UMI="$2"; shift 2 ;;
    --min-snp) MIN_SNP="$2"; shift 2 ;;
    --field) FIELD="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
endone

for var in PLP VCF BARCODES OUTDIR; do
  if [[ -z "${!var}" ]]; then
    echo "[ERROR] --${var,,} is required" >&2
    show_help >&2
    exit 1
  fi
done

for f in "$VCF" "$VCF".tbi "$BARCODES"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing file: $f" >&2
    exit 1
  fi
done

mkdir -p "$OUTDIR"

if [[ -z "$SM_LIST" ]]; then
  SM_LIST="$OUTDIR/samples.tmp"
  docker run --rm -u $(id -u):$(id -g) -v "$HOST_ROOT":"$HOST_ROOT" parisa/genotype:2.6 \
    bcftools query -l "$VCF" > "$SM_LIST"
  CLEAN_SM=1
else
  CLEAN_SM=0
fi

set -x
# shellcheck disable=SC2086
docker run --rm -u $(id -u):$(id -g) -v "$HOST_ROOT":"$HOST_ROOT" \
  --entrypoint /opt/conda/bin/popscle "$IMAGE" demuxlet \
    --plp "$PLP" \
    --vcf "$VCF" \
    --field "$FIELD" \
    --group-list "$BARCODES" \
    --sm-list "$SM_LIST" \
    --doublet-prior "$DOUBLETP" \
    --min-total "$MIN_TOTAL" \
    --min-umi "$MIN_UMI" \
    --min-snp "$MIN_SNP" \
    --out "$OUTDIR/demuxlet"
set +x

if [[ $CLEAN_SM -eq 1 ]]; then rm -f "$SM_LIST"; fi

echo "[INFO] Demuxlet outputs in $OUTDIR"
