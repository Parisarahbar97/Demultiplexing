#!/usr/bin/env bash
#
# mis_prepare.sh
# ---------------
# Combine cleaned per-pool VCFs, drop duplicate donors, merge to a single
# multi-sample dataset and emit chr-prefixed, per-chromosome slices that can be
# uploaded to the Michigan Imputation Server (hg38 panels such as 1000G phase3).
#
# Usage:
#   ./scripts/mis_prepare.sh \
#       --input-dir  imputation_work/01_plink_qc \
#       --output-dir imputation_work/02_mis_prep \
#       --threads    4
#
# Requirements: bcftools (>=1.10).  Run inside Docker/Dockerfile.impute or any
# environment that provides the same toolchain.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: mis_prepare.sh [options]

Required options:
  -i, --input-dir  DIR   Directory containing per-pool *.clean.vcf.gz files.
  -o, --output-dir DIR   Directory for MIS-ready outputs.

Optional:
  -t, --threads    INT   Threads for bcftools merge (default: 4).
  -h, --help             Show this help and exit.

Outputs (under --output-dir):
  unique_vcfs.list              List of per-pool unique VCFs used.
  vcf_unique_manifest.tsv       Mapping of source VCF → retained donor IDs.
  all_donors.samples.txt        Combined donor list.
  all_donors.hg38.chr.vcf.gz    Merged multi-sample VCF (chr1..chr22, VCFv4.2).
  by_chrom/all_donors.chr*.vcf.gz  Per-chromosome files ready for MIS upload.
USAGE
}

input_dir=""
output_dir=""
threads=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-dir)  input_dir="$2"; shift 2 ;;
    -o|--output-dir) output_dir="$2"; shift 2 ;;
    -t|--threads)    threads="$2"; shift 2 ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "${input_dir}"  ]] && { echo "ERROR: --input-dir is required" >&2; usage; exit 1; }
[[ -z "${output_dir}" ]] && { echo "ERROR: --output-dir is required" >&2; usage; exit 1; }

input_dir=$(realpath "${input_dir}")
output_dir=$(realpath -m "${output_dir}")

[[ -d "${input_dir}" ]] || { echo "ERROR: input directory not found: ${input_dir}" >&2; exit 1; }

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found on PATH" >&2; exit 1; }

mkdir -p "${output_dir}"
mkdir -p "${output_dir}/by_chrom"

manifest="${output_dir}/vcf_unique_manifest.tsv"
unique_list="${output_dir}/unique_vcfs.list"
samples_all="${output_dir}/all_donors.samples.txt"
merged_vcf="${output_dir}/all_donors.hg38.raw.vcf.gz"
merged_chr_vcf="${output_dir}/all_donors.hg38.chr.vcf.gz"

rm -f "${manifest}" "${unique_list}" "${samples_all}"
touch "${manifest}" "${samples_all}"

declare -A seen_samples=()

echo "[mis_prepare] Scanning cleaned VCFs in ${input_dir}"
mapfile -t vcfs < <(find "${input_dir}" -maxdepth 1 -type f -name '*.clean.vcf.gz' | sort)
if [[ "${#vcfs[@]}" -eq 0 ]]; then
  echo "ERROR: No *.clean.vcf.gz files found in ${input_dir}" >&2
  exit 1
fi

for vcf in "${vcfs[@]}"; do
  base=$(basename "${vcf}")
  keep_file="${output_dir}/${base}.keep.txt"
  tmp_unique="${output_dir}/${base}.unique.vcf.gz"

  : > "${keep_file}"

  while read -r sm; do
    [[ -z "${sm}" ]] && continue
    if [[ -z "${seen_samples[${sm}]:-}" ]]; then
      seen_samples["${sm}"]=1
      echo "${sm}" >> "${keep_file}"
      printf "%s\t%s\n" "${base}" "${sm}" >> "${manifest}"
      printf "%s\n" "${sm}" >> "${samples_all}"
    fi
  done < <(bcftools query -l "${vcf}")

  if [[ -s "${keep_file}" ]]; then
    bcftools view -S "${keep_file}" "${vcf}" -Oz -o "${tmp_unique}"
    tabix -f -p vcf "${tmp_unique}"
    echo "${tmp_unique}" >> "${unique_list}"
  else
    rm -f "${keep_file}"
  fi
done

sort -u -o "${samples_all}" "${samples_all}"
sample_count=$(wc -l < "${samples_all}")

mapfile -t unique_vcfs < "${unique_list}"
if [[ "${#unique_vcfs[@]}" -eq 0 ]]; then
  echo "ERROR: No unique VCFs generated – check manifest overlap" >&2
  exit 1
fi

echo "[mis_prepare] Retained ${sample_count} unique donor(s) across ${#unique_vcfs[@]} VCFs."

echo "[mis_prepare] Merging with bcftools merge…"
bcftools merge --threads "${threads}" --merge all --missing-to-ref \
  -Oz -o "${merged_vcf}" \
  "${unique_vcfs[@]}"
tabix -f -p vcf "${merged_vcf}"

chr_map="${output_dir}/num_to_chr.map"
seq 1 22 | awk '{printf("%s\tchr%s\n",$1,$1)}' > "${chr_map}"
bcftools annotate --rename-chrs "${chr_map}" "${merged_vcf}" -Oz -o "${merged_chr_vcf}"
tabix -f -p vcf "${merged_chr_vcf}"

# Enforce VCFv4.2 header
hdr_tmp=$(mktemp)
bcftools view -h "${merged_chr_vcf}" > "${hdr_tmp}"
sed -i 's/VCFv4\.3/VCFv4.2/' "${hdr_tmp}"
bcftools reheader -h "${hdr_tmp}" "${merged_chr_vcf}" -o "${merged_chr_vcf}.tmp"
mv -f "${merged_chr_vcf}.tmp" "${merged_chr_vcf}"
tabix -f -p vcf "${merged_chr_vcf}"
rm -f "${hdr_tmp}"

echo "[mis_prepare] Splitting per chromosome…"
for c in $(seq 1 22); do
  out_chr="${output_dir}/by_chrom/all_donors.chr${c}.vcf.gz"
  bcftools view -r "chr${c}" "${merged_chr_vcf}" -Oz -o "${out_chr}"
  tabix -f -p vcf "${out_chr}"
  hdr_tmp=$(mktemp)
  bcftools view -h "${out_chr}" > "${hdr_tmp}"
  sed -i 's/VCFv4\.3/VCFv4.2/' "${hdr_tmp}"
  bcftools reheader -h "${hdr_tmp}" "${out_chr}" -o "${out_chr}.tmp"
  mv -f "${out_chr}.tmp" "${out_chr}"
  tabix -f -p vcf "${out_chr}"
  rm -f "${hdr_tmp}"
done

echo "[mis_prepare] Summary:"
echo "  Donors: ${sample_count}"
printf "  Variants per chromosome:\n"
for f in "${output_dir}"/by_chrom/all_donors.chr*.vcf.gz; do
  chr=$(basename "${f}" .vcf.gz)
  count=$(bcftools view -H "${f}" | wc -l)
  printf "    %s\t%s\n" "${chr}" "${count}"
done

echo "[mis_prepare] Outputs ready in ${output_dir}/by_chrom"
