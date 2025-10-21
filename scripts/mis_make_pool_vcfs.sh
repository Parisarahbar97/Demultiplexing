#!/usr/bin/env bash
#
# mis_make_pool_vcfs.sh
# ----------------------
# Regenerate post-imputation per-pool VCFs from a merged MIS output VCF using
# authoritative donor lists for each library. Ensures every expected donor is
# present; otherwise aborts so the user can restore the missing genotypes
# before demultiplexing.
#
# Usage:
#   ./scripts/mis_make_pool_vcfs.sh \
#       --merged-vcf  /path/to/job1.R2filt.vcf.gz \
#       --lists-dir   /path/to/pool_lists \
#       --output-dir  /path/to/mis_job1_pools_fix
#
# List files must match "*_pool_donors*.txt" by default and contain one donor
# ID per line (additional whitespace and blank/comment lines are ignored).
# Outputs are named "{POOL}.imputed.R2filt.vcf.gz" plus matching ".donors.txt".

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: mis_make_pool_vcfs.sh [options]

Required:
  -m, --merged-vcf  FILE   Merged MIS VCF (e.g. job1.R2filt.vcf.gz).
  -l, --lists-dir   DIR    Directory of pool donor lists (*_pool_donors*.txt).
  -o, --output-dir  DIR    Destination for per-pool VCFs.

Optional:
  -p, --pattern     GLOB   Glob for donor list filenames (default: *_pool_donors*.txt).
  -h, --help               Show this message and exit.

Each donor list file should contain one donor ID per line. Lines beginning
with '#' or blank lines are ignored. The script aborts if any expected donor
is missing from the merged VCF to avoid silently omitting pools from demuxlet.
USAGE
}

merged_vcf=""
lists_dir=""
output_dir=""
pattern="*_pool_donors*.txt"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--merged-vcf) merged_vcf="$2"; shift 2 ;;
    -l|--lists-dir)  lists_dir="$2"; shift 2 ;;
    -o|--output-dir) output_dir="$2"; shift 2 ;;
    -p|--pattern)    pattern="$2"; shift 2 ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "ERROR: Unknown option $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "${merged_vcf}" ]] && { echo "ERROR: --merged-vcf is required" >&2; usage; exit 1; }
[[ -z "${lists_dir}"  ]] && { echo "ERROR: --lists-dir is required" >&2; usage; exit 1; }
[[ -z "${output_dir}" ]] && { echo "ERROR: --output-dir is required" >&2; usage; exit 1; }

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found on PATH" >&2; exit 1; }
command -v tabix >/dev/null 2>&1    || { echo "ERROR: tabix not found on PATH" >&2; exit 1; }

merged_vcf=$(realpath "${merged_vcf}")
lists_dir=$(realpath "${lists_dir}")
output_dir=$(realpath -m "${output_dir}")

[[ -f "${merged_vcf}" ]] || { echo "ERROR: merged VCF not found: ${merged_vcf}" >&2; exit 1; }
[[ -d "${lists_dir}" ]]  || { echo "ERROR: lists directory not found: ${lists_dir}" >&2; exit 1; }

mkdir -p "${output_dir}"

merged_samples=$(mktemp)
tmp_files=("${merged_samples}")
cleanup() { rm -f "${tmp_files[@]}"; }
trap cleanup EXIT

bcftools query -l "${merged_vcf}" | sort > "${merged_samples}"

mapfile -t list_files < <(find "${lists_dir}" -maxdepth 1 -type f -name "${pattern}" | sort)
if [[ "${#list_files[@]}" -eq 0 ]]; then
  echo "ERROR: No donor list files matching '${pattern}' in ${lists_dir}" >&2
  exit 1
fi

echo "[mis_make_pool_vcfs] Using merged VCF: ${merged_vcf}"
echo "[mis_make_pool_vcfs] Writing outputs to: ${output_dir}"

for list_file in "${list_files[@]}"; do
  base=$(basename "${list_file}")
  pool="${base%%_pool_donors*}"

  ordered_tmp=$(mktemp)
  sorted_tmp=$(mktemp)
  missing_tmp=$(mktemp)
  tmp_files+=("${ordered_tmp}" "${sorted_tmp}" "${missing_tmp}")

  # Normalise donor list: strip comments/blank lines and keep original order
  awk 'NF && $1 !~ /^#/ { if (!seen[$1]++) print $1 }' "${list_file}" > "${ordered_tmp}"
  if [[ ! -s "${ordered_tmp}" ]]; then
    echo "[mis_make_pool_vcfs] WARNING: ${base} is empty; skipping ${pool}" >&2
    continue
  fi
  sort "${ordered_tmp}" > "${sorted_tmp}"

  comm -13 "${merged_samples}" "${sorted_tmp}" > "${missing_tmp}"
  if [[ -s "${missing_tmp}" ]]; then
    echo "ERROR: Missing donor(s) for pool ${pool}: $(tr '\n' ' ' < "${missing_tmp}")" >&2
    echo "Ensure all donors were included in MIS prep/postprocess before rerunning." >&2
    exit 1
  fi

  out_vcf="${output_dir}/${pool}.imputed.R2filt.vcf.gz"
  tmp_vcf="${out_vcf}.tmp"
  rm -f "${tmp_vcf}"

  echo "[mis_make_pool_vcfs] Subsetting ${pool} → $(basename "${out_vcf}")"
  bcftools view -S "${ordered_tmp}" -m2 -M2 -v snps "${merged_vcf}" -Oz -o "${tmp_vcf}"
  mv -f "${tmp_vcf}" "${out_vcf}"
  tabix -f -p vcf "${out_vcf}"

  donors_out="${output_dir}/${pool}.donors.txt"
  bcftools query -l "${out_vcf}" > "${donors_out}"
done

count=$(find "${output_dir}" -maxdepth 1 -type f -name "*.imputed.R2filt.vcf.gz" | wc -l | tr -d ' ')
echo "[mis_make_pool_vcfs] Completed ${count} pool VCF(s)."
