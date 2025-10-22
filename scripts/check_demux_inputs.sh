#!/usr/bin/env bash
#
# check_demux_inputs.sh
# ---------------------
# Validate that every sample listed in a Nextflow sample sheet has matching
# resources for demuxlet: BAM, barcodes, imputed VCF, donor list, and that
# contig naming between the BAM and VCF is consistent. Designed to be executed
# inside the Docker images that already provide bcftools and samtools.
#
# Usage:
#   ./scripts/check_demux_inputs.sh [samples.csv] [--path-map HOST:CONTAINER]
#
# Exits non-zero if any sample fails validation and emits a tab-delimited
# summary of the checks performed.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: check_demux_inputs.sh [samples.csv] [--path-map HOST:CONTAINER]

Options:
  samples.csv              CSV with columns sample,bam,barcodes,vcf,sm_list.
  --path-map HOST:CONT     Rewrite path prefix HOST → CONT before accessing files.
  -h, --help               Show this message and exit.

Use --path-map when running inside a container that mounts the host filesystem
at an alternate location (e.g. --path-map /home/pr422:/host).
USAGE
}

sheet="examples/samples.csv"
host_prefix=""
container_prefix=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --path-map)
      [[ $# -ge 2 ]] || { echo "ERROR: --path-map requires HOST:CONTAINER" >&2; exit 1; }
      host_prefix=${2%%:*}
      container_prefix=${2#*:}
      shift 2
      ;;
    *)
      sheet="$1"
      shift
      ;;
  esac
done

map_path() {
  local path="$1"
  if [[ -n "${host_prefix}" && -n "${container_prefix}" && "${path}" == "${host_prefix}"* ]]; then
    printf "%s%s" "${container_prefix}" "${path#${host_prefix}}"
  else
    printf "%s" "${path}"
  fi
}

[[ -f "${sheet}" ]] || { echo "ERROR: sample sheet not found: ${sheet}" >&2; exit 1; }

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not on PATH" >&2; exit 1; }
command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not on PATH"  >&2; exit 1; }

tmpdir=$(mktemp -d)
cleanup() { rm -rf "${tmpdir}"; }
trap cleanup EXIT

printf "sample\tstatus\tmessage\n"

failures=0

tail -n +2 "${sheet}" > "${tmpdir}/rows.csv"
while IFS=, read -r sample bam barcodes vcf sm_list; do
  sample=${sample//$'\r'/}
  bam=${bam//$'\r'/}
  barcodes=${barcodes//$'\r'/}
  vcf=${vcf//$'\r'/}
  sm_list=${sm_list//$'\r'/}

  if [[ -z "${sample}" ]]; then
    continue
  fi

  messages=()
  status="OK"

  bam_path=$(map_path "${bam}")
  barcodes_path=$(map_path "${barcodes}")
  vcf_path=$(map_path "${vcf}")
  sm_list_path=$(map_path "${sm_list}")

  if [[ ! -f "${bam_path}" ]]; then
    messages+=("missing BAM")
  fi
  if [[ ! -f "${barcodes_path}" ]]; then
    messages+=("missing barcodes")
  fi
  if [[ ! -f "${vcf_path}" ]]; then
    messages+=("missing VCF")
  fi
  if [[ -z "${sm_list}" || ! -f "${sm_list_path}" ]]; then
    messages+=("missing donor list")
  fi

  if [[ "${#messages[@]}" -eq 0 ]]; then
    donors_expected="${tmpdir}/${sample}.donors.expected"
    donors_found="${tmpdir}/${sample}.donors.found"
    comm_missing="${tmpdir}/${sample}.donors.missing"
    comm_extra="${tmpdir}/${sample}.donors.extra"

    awk 'NF && $1 !~ /^#/' "${sm_list_path}" | sort -u > "${donors_expected}"
    bcftools query -l "${vcf_path}" | sort -u > "${donors_found}"

    comm -23 "${donors_expected}" "${donors_found}" > "${comm_missing}"
    comm -13 "${donors_expected}" "${donors_found}" > "${comm_extra}"

    if [[ -s "${comm_missing}" ]]; then
      messages+=("donors missing from VCF: $(tr '\n' ' ' < "${comm_missing}")")
    fi
    if [[ -s "${comm_extra}" ]]; then
      messages+=("unexpected donors in VCF: $(tr '\n' ' ' < "${comm_extra}")")
    fi

    bam_contigs="${tmpdir}/${sample}.bam.contigs"
    vcf_contigs="${tmpdir}/${sample}.vcf.contigs"
    diff_missing="${tmpdir}/${sample}.contigs.missing"
    diff_extra="${tmpdir}/${sample}.contigs.extra"

    samtools idxstats "${bam_path}" \
      | awk '$1 != "*" && $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/{print $1}' \
      | sort -u > "${bam_contigs}"
    bcftools view -h "${vcf_path}" \
      | awk -F'[=,]' '/^##contig=/{gsub(/ID=/,"",$2); if ($2 ~ /^chr([1-9]|1[0-9]|2[0-2])$/) print $2}' \
      | sort -u > "${vcf_contigs}"

    comm -23 "${bam_contigs}" "${vcf_contigs}" > "${diff_missing}"
    comm -13 "${bam_contigs}" "${vcf_contigs}" > "${diff_extra}"

    if [[ -s "${diff_missing}" ]]; then
      messages+=("VCF missing contigs: $(tr '\n' ' ' < "${diff_missing}")")
    fi
    if [[ -s "${diff_extra}" ]]; then
      messages+=("VCF has extra contigs: $(tr '\n' ' ' < "${diff_extra}")")
    fi
  fi

  if [[ "${#messages[@]}" -gt 0 ]]; then
    status="FAIL"
    printf "%s\t%s\t%s\n" "${sample}" "${status}" "$(IFS='; '; echo "${messages[*]}")"
    failures=$((failures + 1))
  else
    printf "%s\t%s\t%s\n" "${sample}" "${status}" "All checks passed"
  fi
done < "${tmpdir}/rows.csv"

if [[ "${failures}" -gt 0 ]]; then
  exit 1
fi

echo "All samples passed validation."
