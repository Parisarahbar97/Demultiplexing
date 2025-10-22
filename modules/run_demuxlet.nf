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
  def prefix = plp_plp.simpleName.replaceFirst(/\.plp$/, '')
  """
  set -euo pipefail
  PREFIX="${prefix}"

  BC_TMP=\$(mktemp)
  # Decompress if needed
  case "$barcodes" in
    *.gz)  gzip -cd "$barcodes" > "\$BC_TMP.raw" ;;
    *)     cat  "$barcodes" > "\$BC_TMP.raw" ;;
  esac
  # If it looks like CSV/TSV with a header, strip it and keep first column
  first_line=\$(head -n1 "\$BC_TMP.raw")
  if echo "\$first_line" | grep -qi 'barcode'; then
    awk -F'[,\t]' 'NR>1{print \$1}' "\$BC_TMP.raw" > "\$BC_TMP"
  else
    mv "\$BC_TMP.raw" "\$BC_TMP"
  fi

  popscle demuxlet \\
    --plp "\$PREFIX" \\
    --vcf "$vcf_for_demux" \\
    --field "${params.field}" \\
    --group-list "\$BC_TMP" \\
    --sm-list "$sm_list" \\
    --doublet-prior ${params.doublet_prior} \\
    --out demuxlet

  rm -f "\$BC_TMP" "\$BC_TMP.raw"
  """
}
