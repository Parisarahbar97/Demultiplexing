process DEMUXLET {
  tag "${sample_id}"
  label 'DEMUXLET'

  input:
  tuple val(sample_id),
       path(plp_plp), path(plp_var), path(plp_cel),
       path(vcf_clean), path(sm_list), path(barcodes)

  output:
  tuple val(sample_id),
       path("${sample_id}_demuxlet/demuxlet.best"),
       path(barcodes), path(sm_list)

  script:
  """
  set -euo pipefail

  # reconstruct prefix from .plp.gz path
  PLP_PREFIX=\$(echo "${plp_plp}" | sed 's/\\.plp\\.gz\$//')

  mkdir -p ${sample_id}_demuxlet

  # demuxlet prefers a plain barcode list
  BC_TXT=\$(mktemp)
  if [[ "${barcodes}" == *.gz ]]; then
    gzip -cd "${barcodes}" > "\$BC_TXT"
  else
    cp "${barcodes}" "\$BC_TXT"
  fi

  popscle demuxlet \\
    --plp "\$PLP_PREFIX" \\
    --vcf "${vcf_clean}" \\
    --field ${params.field} \\
    --group-list "\$BC_TXT" \\
    --sm-list "${sm_list}" \\
    --doublet-prior ${params.doublet_prior} \\
    --out "${sample_id}_demuxlet/demuxlet"
  """
}
