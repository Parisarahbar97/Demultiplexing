#!/usr/bin/env bash
set -euo pipefail

# Usage: prep_healthy_vcfs.sh /in /out [--chr-map /out/chr_map.txt] [--assume-grch38]
# - /in  : folder with your healthy donor VCFs (e.g., D1A.vcf, D1P.vcf, ...)
# - /out : folder to write outputs (will be created if missing)
# - --chr-map: optional 2-column file to rename contigs (e.g., "1 chr1" lines)
# - --assume-grch38: skip REF/ALT normalization against a FASTA (keeps script minimal)

IN="${1:-/in}"
OUT="${2:-/out}"
CHR_MAP=""
ASSUME_GRCH38="no"

shift 2 || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chr-map)        CHR_MAP="$2"; shift 2;;
    --assume-grch38)  ASSUME_GRCH38="yes"; shift 1;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

mkdir -p "$OUT"
echo "Input dir : $IN"
echo "Output dir: $OUT"
echo "chr-map    : ${CHR_MAP:-<none>}"
echo

echo "Step 0) Quick header/contig check on the first VCF..."
FIRST_VCF="$(ls "$IN"/D*.vcf* 2>/dev/null | head -n1 || true)"
if [[ -z "${FIRST_VCF}" ]]; then
  echo "No VCFs like D*.vcf found in $IN"; exit 1
fi
bcftools view -h "$FIRST_VCF" | sed -n '1,30p' || true
echo
echo "Tip: Make sure contig style (e.g., 'chr1' vs '1') matches your BAMs."

echo "Step 1) Sort + bgzip + tabix each donor VCF..."
shopt -s nullglob
mapfile -t VCF_LIST < <(ls "$IN"/D*.vcf "$IN"/D*.vcf.gz 2>/dev/null || true)
if [[ ${#VCF_LIST[@]} -eq 0 ]]; then
  echo "Found 0 donor VCFs matching D*.vcf*"; exit 1
fi

TMP_LIST=()
for V in "${VCF_LIST[@]}"; do
  BASE=$(basename "$V")
  BASE_NOEXT="${BASE%.vcf.gz}"; BASE_NOEXT="${BASE_NOEXT%.vcf}"
  OUTVCF="$OUT/${BASE_NOEXT}.sorted.vcf.gz"
  echo "  - $BASE → $OUTVCF"
  bcftools sort -Oz -o "$OUTVCF" "$V"
  tabix -p vcf "$OUTVCF"
  TMP_LIST+=("$OUTVCF")
done

echo "Step 2) (Optional) Rename chromosomes if contigs don't match BAMs..."
if [[ -n "$CHR_MAP" ]]; then
  echo "  Applying --chr-map $CHR_MAP to each file..."
  RENAMED_LIST=()
  for V in "${TMP_LIST[@]}"; do
    OUTVCF="${V%.vcf.gz}.renamed.vcf.gz"
    bcftools annotate --rename-chrs "$CHR_MAP" "$V" -Oz -o "$OUTVCF"
    tabix -p vcf "$OUTVCF"
    RENAMED_LIST+=("$OUTVCF")
  done
  TMP_LIST=("${RENAMED_LIST[@]}")
else
  echo "  Skipping rename (no --chr-map provided)."
fi

echo "Step 3) Merge all donors → one multi-sample VCF..."
MERGED="$OUT/healthy_all.merged.vcf.gz"
bcftools merge -Oz -o "$MERGED" "${TMP_LIST[@]}"
tabix -p vcf "$MERGED"

echo "Step 4) Keep only biallelic SNPs (cleanest for demux)..."
BIALLELIC="$OUT/healthy_all.merged.biallelic.vcf.gz"
bcftools view -m2 -M2 -v snps -Oz -o "$BIALLELIC" "$MERGED"
tabix -p vcf "$BIALLELIC"

echo "Step 5) Add AF, AC, AN tags to INFO (handy for checks/pileups)..."
FOR_PILEUP="$OUT/healthy_all.for_pileup.vcf.gz"
bcftools +fill-tags "$BIALLELIC" -- -t AC,AN,AF -Oz -o "$FOR_PILEUP"
tabix -p vcf "$FOR_PILEUP"

echo "Step 6) Write out sample list (donor IDs in the merged VCF header)..."
bcftools query -l "$BIALLELIC" > "$OUT/healthy_all.samples.txt"

echo
echo "Done."
echo "Outputs:"
echo "  - $BIALLELIC            (multi-sample, biallelic SNPs)  ← use for demuxlet"
echo "  - $FOR_PILEUP           (same sites + AF/AC/AN)         ← use for dsc-pileup"
echo "  - $OUT/healthy_all.samples.txt  (check donor IDs)"
