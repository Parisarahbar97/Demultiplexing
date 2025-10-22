process dsc_pileup {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'

  input:
  tuple val(sample_id), path(bam), path(barcodes), path(vcf_for_pileup)

  output:
  tuple val(sample_id),
        path("${sample_id}_pileup.plp.gz"),
        path("${sample_id}_pileup.var.gz"),
        path("${sample_id}_pileup.cel.gz")

  script:
  def pileup_opts = []
    if( params.min_mapq != null )  pileup_opts << "--min-MQ ${params.min_mapq}"
    if( params.min_baseq != null ) pileup_opts << "--min-BQ ${params.min_baseq}"
    if( params.cap_bq   != null )  pileup_opts << "--cap-BQ ${params.cap_bq}"
  def pileup_opts_str = pileup_opts.join(' ')
  """
  set -euo pipefail
 BC_TMP=$(mktemp)
# Decompress if needed
case "$barcodes" in
  *.gz)  zcat "$barcodes" > "$BC_TMP.raw" ;;
  *)     cat  "$barcodes" > "$BC_TMP.raw" ;;
esac
# If it looks like CSV with a header, strip header and first column only
first_line=$(head -n1 "$BC_TMP.raw")
if echo "$first_line" | grep -qi 'barcode'; then
  awk -F'[,\t]' 'NR>1{print $1}' "$BC_TMP.raw" > "$BC_TMP"
else
  # Already 1 barcode per line
  mv "$BC_TMP.raw" "$BC_TMP"
fi

  popscle dsc-pileup \\
    --sam "$bam" \\
    --vcf "$vcf_for_pileup" \\
    --group-list "\$BC_TMP" \\
    --tag-group "${params.tag_group}" \\
    --tag-UMI "${params.tag_umi}" \\
    --out "${sample_id}_pileup" ${pileup_opts_str}

  rm -f "\$BC_TMP"
  """
}
