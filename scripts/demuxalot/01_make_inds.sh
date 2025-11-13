#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 01_make_inds.sh --vcf INPUT --out OUTPUT [options]

Create a Demuxalot sample list (one donor ID per line) by querying the VCF.
The VCF must be bgzip-compressed with a matching index.

Options:
  --vcf PATH          Input VCF/BCF (bgzip + index)
  --out PATH          Output text file (one sample ID per line)
  --host-root PATH    Host directory to bind inside Docker (default: /home/pr422)
  --image NAME        Docker image containing bcftools (default: parisa/genotype:2.6)
  -h, --help          Show this help message
USAGE
}

VCF=""
OUT=""
HOST_ROOT=${HOST_ROOT:-/home/pr422}
IMAGE=${GENO_IMAGE:-parisa/genotype:2.6}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

if [[ -z "$VCF" || -z "$OUT" ]]; then
  echo "[ERROR] --vcf and --out are required" >&2
  show_help >&2
  exit 1
fi

if [[ ! -f "$VCF" ]]; then
  echo "[ERROR] Missing VCF: $VCF" >&2
  exit 1
fi

if [[ ! -f "$VCF.tbi" && ! -f "$VCF.csi" ]]; then
  echo "[ERROR] Missing VCF index (.tbi or .csi) for $VCF" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT")"

# shellcheck disable=SC2086
docker run --rm -u "$(id -u)":"$(id -g)" \
  -v "$HOST_ROOT":"$HOST_ROOT" \
  "$IMAGE" bcftools query -l "$VCF" > "$OUT"

if [[ ! -s "$OUT" ]]; then
  echo "[WARN] $OUT is empty; did bcftools query detect any samples?" >&2
else
  echo "[INFO] Wrote $(wc -l < "$OUT") sample IDs to $OUT"
fi
