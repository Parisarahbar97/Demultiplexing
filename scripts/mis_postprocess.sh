#!/usr/bin/env bash
#
# mis_postprocess.sh
# -------------------
# Post-process Michigan Imputation Server results by converting each chr_*.zip
# bundle to BCF, concatenating across chromosomes, filtering by INFO/R2, and
# optionally applying a minor allele frequency window. Outputs a final VCF
# suitable for demuxlet (GP/DS preserved).
#
# Usage:
#   ./scripts/mis_postprocess.sh \
#       --input-dir  imputation_work/03_imputed/mis_job1_raw \
#       --output-dir imputation_work/03_imputed/mis_job1_post \
#       --prefix     job1 \
#       --r2-min     0.4 \
#       --maf-min    0.0 \
#       --maf-max    1.0
#
# Requirements: unzip, bcftools (>=1.10), tabix.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: mis_postprocess.sh [options]

Required:
  -i, --input-dir  DIR    Directory containing chr_*.zip from MIS.
  -o, --output-dir DIR    Directory for processed outputs.

Optional:
  -p, --prefix     STR    Prefix for output files (default: imputed).
  -r, --r2-min     FLOAT  Minimum INFO/R2 to retain (default: 0.4).
  -m, --maf-min    FLOAT  Minimum AF to retain (default: 0.0).
  -M, --maf-max    FLOAT  Maximum AF to retain (default: 1.0).
  -h, --help              Show this help and exit.

Outputs (under --output-dir):
  bcf/{chrN.bcf}                    Intermediate BCFs per chromosome.
  post/{prefix}.merged.bcf          Concatenated BCF across chr1..chr22.
  post/{prefix}.R2filt.vcf.gz       After INFO/R2 filter (autosomal SNPs).
  post/{prefix}.R2filt.MAF.vcf.gz   After INFO/R2 + AF window.
USAGE
}

input_dir=""
output_dir=""
prefix="imputed"
r2_min="0.4"
maf_min="0.0"
maf_max="1.0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-dir)  input_dir="$2"; shift 2 ;;
    -o|--output-dir) output_dir="$2"; shift 2 ;;
    -p|--prefix)     prefix="$2"; shift 2 ;;
    -r|--r2-min)     r2_min="$2"; shift 2 ;;
    -m|--maf-min)    maf_min="$2"; shift 2 ;;
    -M|--maf-max)    maf_max="$2"; shift 2 ;;
        --zip-password) zip_password="$2"; shift 2 ;;
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

unz_dir="${output_dir}/unzipped"
bcf_dir="${output_dir}/bcf"
post_dir="${output_dir}/post"
mkdir -p "${unz_dir}" "${bcf_dir}" "${post_dir}"
rm -f "${bcf_dir}"/*.bcf "${bcf_dir}"/*.bcf.csi 2>/dev/null || true

echo "[mis_postprocess] Converting chr_*.zip to BCF"
for zip_file in "${input_dir}"/chr_*.zip; do
  [[ -f "${zip_file}" ]] || continue
  if [[ -n "${zip_password}" ]]; then
    unzip -P "${zip_password}" -o "${zip_file}" -d "${unz_dir}" >/dev/null
  else
    unzip -o "${zip_file}" -d "${unz_dir}" >/dev/null
  fi
  base=$(basename "${zip_file}" .zip)
  chr=${base#chr_}
  src="${unz_dir}/chr${chr}.dose.vcf.gz"
  if [[ ! -f "${src}" ]]; then
    echo "ERROR: expected file ${src} is missing" >&2
    exit 1
  fi
  bcftools view -Ob -o "${bcf_dir}/chr${chr}.bcf" "${src}"
  bcftools index -f "${bcf_dir}/chr${chr}.bcf"
  rm -f "${src}" "${unz_dir}/chr${chr}.info.gz"
done

echo "[mis_postprocess] Concatenating chr1..chr22"
merged_bcf="${post_dir}/${prefix}.merged.bcf"
bcftools concat -Ob -o "${merged_bcf}" \
  "${bcf_dir}"/chr1.bcf "${bcf_dir}"/chr2.bcf "${bcf_dir}"/chr3.bcf "${bcf_dir}"/chr4.bcf \
  "${bcf_dir}"/chr5.bcf "${bcf_dir}"/chr6.bcf "${bcf_dir}"/chr7.bcf "${bcf_dir}"/chr8.bcf \
  "${bcf_dir}"/chr9.bcf "${bcf_dir}"/chr10.bcf "${bcf_dir}"/chr11.bcf "${bcf_dir}"/chr12.bcf \
  "${bcf_dir}"/chr13.bcf "${bcf_dir}"/chr14.bcf "${bcf_dir}"/chr15.bcf "${bcf_dir}"/chr16.bcf \
  "${bcf_dir}"/chr17.bcf "${bcf_dir}"/chr18.bcf "${bcf_dir}"/chr19.bcf "${bcf_dir}"/chr20.bcf \
  "${bcf_dir}"/chr21.bcf "${bcf_dir}"/chr22.bcf
bcftools index -f "${merged_bcf}"

echo "[mis_postprocess] Filtering INFO/R2 >= ${r2_min}"
r2_vcf="${post_dir}/${prefix}.R2filt.vcf.gz"
bcftools view -i "INFO/R2>=${r2_min}" -v snps -m2 -M2 "${merged_bcf}" -Ou \
  | bcftools sort -Oz -o "${r2_vcf}"
tabix -f -p vcf "${r2_vcf}"

echo "[mis_postprocess] Applying AF window ${maf_min}-${maf_max}"
maf_vcf="${post_dir}/${prefix}.R2filt.MAF.vcf.gz"
bcftools +fill-tags "${r2_vcf}" -Ou -- -t AF \
  | bcftools view -i "AF>=${maf_min} && AF<=${maf_max}" -Oz -o "${maf_vcf}"
tabix -f -p vcf "${maf_vcf}"

printf "[mis_postprocess] Summary:\n"
printf "  Samples: %s\n" "$(bcftools query -l "${r2_vcf}" | wc -l)"
printf "  Variants after R2: %s\n" "$(bcftools view -H "${r2_vcf}" | wc -l)"
printf "  Variants after R2+MAF: %s\n" "$(bcftools view -H "${maf_vcf}" | wc -l)"

echo "[mis_postprocess] Final VCF ready for demuxlet: ${maf_vcf}"
