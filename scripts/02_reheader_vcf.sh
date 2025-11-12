#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 02_reheader_vcf.sh --vcf INPUT --bam INPUT --out OUTPUT [options]

Reheader a pool VCF so its contig order matches the BAM header. Uses bcftools
+ samtools inside Docker (defaults defined below).

Options:
  --vcf PATH          Input VCF/BCF (bgzip + index required)
  --bam PATH          Coordinate-sorted BAM (with .bai)
  --out PATH          Output VCF path (bgzip + tabix index created)
  --host-root PATH    Host path to mount inside Docker (default: /home/pr422)
  --image NAME        Docker image with bcftools/samtools (default: parisa/genotype:2.6)
  -h, --help          Show usage
USAGE
}

VCF=""
BAM=""
OUT=""
HOST_ROOT=${HOST_ROOT:-/home/pr422}
IMAGE=${GENO_IMAGE:-parisa/genotype:2.6}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2 ;;
    --bam) BAM="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

if [[ -z "$VCF" || -z "$BAM" || -z "$OUT" ]]; then
  echo "[ERROR] --vcf, --bam, and --out are required" >&2
  show_help >&2
  exit 1
fi

if [[ ! -f "$VCF" ]]; then
  echo "[ERROR] Missing VCF: $VCF" >&2
  exit 1
fi

if [[ -f "$VCF.tbi" ]]; then
  VCF_IDX="$VCF.tbi"
elif [[ -f "$VCF.csi" ]]; then
  VCF_IDX="$VCF.csi"
else
  echo "[ERROR] Missing VCF index (.tbi or .csi) for $VCF" >&2
  exit 1
fi

if [[ ! -f "$BAM" || ! -f "$BAM.bai" ]]; then
  echo "[ERROR] Missing BAM or BAM index for $BAM" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT")"

export VCF BAM OUT VCF_IDX

read -r -d '' CMD <<'SCRIPT' || true
set -euo pipefail
IN_VCF="$VCF"
IN_VCF_IDX="$VCF_IDX"
IN_BAM="$BAM"
OUT_VCF="$OUT"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

> "$TMP_DIR/contigs.tsv"
while IFS=$'\t' read -r -a fields; do
  label=${fields[0]:-}
  if [[ $label == "@SQ" ]]; then
    contig=""; length=""
    for field in "${fields[@]:1}"; do
      case $field in
        SN:*) contig=${field#SN:} ;;
        LN:*) length=${field#LN:} ;;
      esac
    done
    if [[ -n $contig ]]; then
      [[ -z $length ]] && length=1
      printf '%s\t%s\n' "$contig" "$length" >> "$TMP_DIR/contigs.tsv"
    fi
  fi
done < <(samtools view -H "$IN_BAM")

awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' "$TMP_DIR/contigs.tsv" > "$TMP_DIR/contigs_header.txt"

bcftools view -h "$IN_VCF" > "$TMP_DIR/orig_header.txt"

awk 'NR==FNR {contig_lines[++n]=$0; next} /^##contig=/ {next} /^#CHROM/ {for(i=1;i<=n;i++) print contig_lines[i]; print $0; next} {print}' \
  "$TMP_DIR/contigs_header.txt" "$TMP_DIR/orig_header.txt" > "$TMP_DIR/new_header.txt"

bcftools index -f "$IN_VCF"
cp "$IN_VCF" "$TMP_DIR/input.vcf.gz"
cp "$IN_VCF_IDX" "$TMP_DIR/input.vcf.gz.$(basename "$IN_VCF_IDX" | sed 's/.*\.//')" 2>/dev/null || true

bcftools reheader -h "$TMP_DIR/new_header.txt" "$TMP_DIR/input.vcf.gz" -o "$TMP_DIR/reheader.vcf.gz"
tabix -f -p vcf "$TMP_DIR/reheader.vcf.gz"

bcftools sort -O z -o "$OUT_VCF" "$TMP_DIR/reheader.vcf.gz"
tabix -f -p vcf "$OUT_VCF"
SCRIPT

# shellcheck disable=SC2086
docker run --rm \
  -u $(id -u):$(id -g) \
  -e VCF="$VCF" \
  -e VCF_IDX="$VCF_IDX" \
  -e BAM="$BAM" \
  -e OUT="$OUT" \
  -v "$HOST_ROOT":"$HOST_ROOT" \
  "$IMAGE" bash -lc "$CMD"

echo "[INFO] Wrote $OUT"
