#!/usr/bin/env bash
#
# mis_postprocess.sh
# -------------------
# Post-process Michigan Imputation Server results:
#   1. Unzip per-chromosome archives.
#   2. Concatenate chr1..chr22 into a single VCF.
#   3. Filter by INFO/R2 and optional minor allele frequency.
#   4. Produce a final VCF ready for demuxlet (GP/DS in FORMAT).
#
# Usage:
#   ./scripts/mis_postprocess.sh \
#       --input-dir  imputation_work/03_imputed/mis_job1_raw \
#       --output-dir imputation_work/03_imputed/mis_job1_post \
#       --prefix     job1 \
#       --r2-min     0.4 \
#       --maf-min    0.05
#
# Requirements: unzip, bcftools, tabix.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: mis_postprocess.sh [options]

Required:
  -i, --input-dir  DIR    Directory containing chr_*.zip from MIS.
  -o, --output-dir DIR    Directory for processed outputs.

Optional:
  -p, --prefix     STR    Prefix for output VCFs (default: imputed).
  -r, --r2-min     FLOAT  Minimum INFO/R2 to retain (default: 0.4).
  -m, --maf-min    FLOAT  Minimum AF to retain (default: 0.05).
  -M, --maf-max    FLOAT  Maximum AF to retain (default: 0.95).
  -h, --help              Show this help and exit.

Outputs:
  merge/{prefix}.merged.vcf.gz            Concatenated chr1..22 VCF.
  post/{prefix}.R2filt.vcf.gz             After INFO/R2 filter.
  post/{prefix}.R2filt.MAF.vcf.gz         After R2 + MAF filters.
USAGE
}

input_dir=""
output_dir=""
prefix="imputed"
r2_min="0.4"
maf_min="0.05"
maf_max="0.95"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-dir)  input_dir="$2"; shift 2 ;;
    -o|--output-dir) output_dir="$2"; shift 2 ;;
    -p|--prefix)     prefix="$2"; shift 2 ;;
    -r|--r2-min)     r2_min="$2"; shift 2 ;;
    -m|--maf-min)    maf_min="$2"; shift 2 ;;
    -M|--maf-max)    maf_max="$2"; shift 2 ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "${input_dir}"  ]] && { echo "ERROR: --input-dir is required" >&2; usage; exit 1; }
[[ -z "${output_dir}" ]] && { echo "ERROR: --output-dir is required" >&2; usage; exit 1; }

input_dir=$(realpath "${input_dir}")
output_dir=$(realpath -m "${output_dir}")

[[ -d "${input_dir}" ]] || { echo "ERROR: input dir not found: ${input_dir}" >&2; exit 1; }

command -v unzip >/dev/null 2>&1     || { echo "ERROR: unzip not found" >&2; exit 1; }
command -v bcftools >/dev/null 2>&1  || { echo "ERROR: bcftools not found" >&2; exit 1; }
command -v tabix >/dev/null 2>&1     || { echo "ERROR: tabix not found" >&2; exit 1; }

mkdir -p "${output_dir}"
unz_dir="${output_dir}/unzip"
merge_dir="${output_dir}/merge"
post_dir="${output_dir}/post"
mkdir -p "${unz_dir}" "${merge_dir}" "${post_dir}"

echo "[mis_postprocess] Unzipping chr_*.zip from ${input_dir}"
for zip_file in "${input_dir}"/chr_*.zip; do
  [[ -f "${zip_file}" ]] || continue
  unzip -o "${zip_file}" -d "${unz_dir}" >/dev/null
done

find "${unz_dir}" -name "*.vcf.gz" | sort > "${merge_dir}/vcfs.list"
if [[ ! -s "${merge_dir}/vcfs.list" ]]; then
  echo "ERROR: No VCF files found after unzip. Check input directory." >&2
  exit 1
fi

echo "[mis_postprocess] Concatenating chr1..chr22"
merged_vcf="${merge_dir}/${prefix}.merged.vcf.gz"
: > "${merge_dir}/chr_vcfs.list"
for c in $(seq 1 22); do
  v=$(grep -m1 -E "/chr[_-]?${c}[^/]*\.vcf\.gz$" "${merge_dir}/vcfs.list" || true)
  if [[ -z "${v}" ]]; then
    echo "ERROR: Missing chromosome ${c} VCF in ${merge_dir}/vcfs.list" >&2
    exit 1
  fi
  echo "${v}" >> "${merge_dir}/chr_vcfs.list"
done
bcftools concat -f "${merge_dir}/chr_vcfs.list" -Oz -o "${merged_vcf}"
tabix -f -p vcf "${merged_vcf}"

echo "[mis_postprocess] Filtering INFO/R2 >= ${r2_min}"
r2_vcf="${post_dir}/${prefix}.R2filt.vcf.gz"
bcftools view -i "INFO/R2>=${r2_min}" -v snps -m2 -M2 -r chr1-chr22 "${merged_vcf}" -Ou \
  | bcftools sort -Oz -o "${r2_vcf}"
tabix -f -p vcf "${r2_vcf}"

echo "[mis_postprocess] Applying AF window ${maf_min}-${maf_max}"
maf_vcf="${post_dir}/${prefix}.R2filt.MAF.vcf.gz"
bcftools +fill-tags "${r2_vcf}" -Ou -- -t AF \
  | bcftools view -i "AF>=${maf_min} && AF<=${maf_max}" -Oz -o "${maf_vcf}"
tabix -f -p vcf "${maf_vcf}"

printf "[mis_postprocess] Summary:\n"
printf "  Samples: %s\n" "$(bcftools query -l "${maf_vcf}" | wc -l)"
printf "  Variants after R2: %s\n" "$(bcftools view -H "${r2_vcf}" | wc -l)"
printf "  Variants after R2+MAF: %s\n" "$(bcftools view -H "${maf_vcf}" | wc -l)"

echo "[mis_postprocess] Final VCF ready for demuxlet: ${maf_vcf}"
