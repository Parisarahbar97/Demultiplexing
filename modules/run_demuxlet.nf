process run_demuxlet {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/${sample_id}_demuxlet", mode: 'copy'

  input:
  tuple val(sample_id),
        path(plp_plp), path(plp_var), path(plp_cel),
        path(vcf_for_demux), path(barcodes), path(sm_list)

  output:
  tuple val(sample_id), path("demuxlet.best")

  script:
  """
  set -euo pipefail
  PREFIX="${sample_id}_pileup"

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

  popscle demuxlet \\
    --plp "\$PREFIX" \\
    --vcf "$vcf_for_demux" \\
    --field "${params.field}" \\
    --group-list "\$BC_TMP" \\
    --sm-list "$sm_list" \\
    --doublet-prior ${params.doublet_prior} \\
    --out demuxlet

  rm -f "\$BC_TMP"
  """
}
