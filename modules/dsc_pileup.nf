process PILEUP {
  tag "${sample_id}"
  label 'PILEUP'

  input:
  tuple val(sample_id), path(bam), path(barcodes), path(vcf_clean), path(sm_list)

  output:
  // pass along everything needed for the next step
  tuple val(sample_id),
       path("${sample_id}_pileup.plp.gz"),
       path("${sample_id}_pileup.var.gz"),
       path("${sample_id}_pileup.cel.gz"),
       path(vcf_clean),
       path(sm_list),
       path(barcodes)

  script:
  """
  set -euo pipefail

  # demuxlet wants a plain list for group-list; unzip if needed
  BC_TXT=\$(mktemp)
  if [[ "${barcodes}" == *.gz ]]; then
    gzip -cd "${barcodes}" > "\$BC_TXT"
  else
    cp "${barcodes}" "\$BC_TXT"
  fi

  popscle dsc-pileup \\
    --sam "${bam}" \\
    --vcf "${vcf_clean}" \\
    --group-list "\$BC_TXT" \\
    --tag-group ${params.tag_group} \\
    --tag-UMI ${params.tag_umi} \\
    --out "${sample_id}_pileup"

  cp "\$BC_TXT" barcodes.txt || true
  """
}
