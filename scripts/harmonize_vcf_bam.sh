#!/usr/bin/env bash
set -euo pipefail

check_exit_codes () {
  local GET_PIPESTATUS="${PIPESTATUS[@]}"; local ec
  for ec in ${GET_PIPESTATUS}; do
    if [ "${ec}" -ne 0 ]; then return "${ec}"; fi
  done
  return 0
}

check_if_programs_exists () {
  command -v awk >/dev/null 2>&1 || { echo 'Error: "awk" not found.' >&2; return 2; }
  command -v bcftools >/dev/null 2>&1 || { echo 'Error: "bcftools" not found.' >&2; return 2; }
  command -v samtools >/dev/null 2>&1 || { echo 'Error: "samtools" not found.' >&2; return 2; }
  return 0
}

get_contig_order_from_bam () {
  local bam_input_file="${1}"
  local output_type="${2}"

  if [ $# -ne 2 ]; then
    printf 'Usage: get_contig_order_from_bam BAM_file output_type (names|chrom_sizes|vcf)\n' >&2
    return 1
  fi

  check_if_programs_exists || return $?

  samtools view -H "${bam_input_file}" \
  | awk -F '\t' -v out="${output_type}" '
    $1=="@SQ"{
      n+=1; name=""; len="";
      for(i=2;i<=NF;i++){
        if ($i ~ /^SN:/){ name=substr($i,4) }
        else if ($i ~ /^LN:/){ len=substr($i,4) }
      }
      nm[n]=name; ln[n]=len
    }
    END{
      if (n==0){ print "Error: No @SQ in BAM header." > "/dev/stderr"; exit 1 }
      if (out=="names"){
        for(i=1;i<=n;i++){ printf("%s%s", (i>1?" ":""), nm[i]) } print ""
      } else if (out=="chrom_sizes"){
        for(i=1;i<=n;i++){ print nm[i] "\t" ln[i] }
      } else if (out=="vcf"){
        for(i=1;i<=n;i++){ print "##contig=<ID=" nm[i] ",length=" ln[i] ">" }
      } else {
        print "Error: unknown output_type " out > "/dev/stderr"; exit 1
      }
    }'
  check_exit_codes; return $?
}

sort_vcf_same_as_bam () {
  local bam_input_file="${1}"
  local vcf_input_file="${2}"
  local vcf_type="${3:-v}"   # v|z|u|b

  if [ $# -lt 2 ]; then
    printf 'Usage: sort_vcf_same_as_bam BAM_file VCF_file [v|z|u|b]\n' >&2
    return 1
  fi

  check_if_programs_exists || return $?

  cat \
    <( bcftools view -h "${vcf_input_file}" \
       | awk '{ if ($0 !~ /^##contig=/ && $0 !~ /^#CHROM/) print }' \
       ; get_contig_order_from_bam "${bam_input_file}" vcf \
       ; bcftools view -h "${vcf_input_file}" | tail -n 1 \
     ) \
    <( bcftools view -H -O v "${vcf_input_file}" ) \
  | bcftools sort -O "${vcf_type}"
  check_exit_codes; return $?
}

sort_vcf_same_as_bam "$@"
