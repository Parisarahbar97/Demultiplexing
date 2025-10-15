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
  PREFIX=\$(echo "${plp_plp.simpleName}" | sed 's/\\.plp\$//')

  BC_TMP=\$(mktemp)
  case "$barcodes" in
    *.gz) gzip -cd "$barcodes" > "\$BC_TMP" ;;
    *)    cp "$barcodes" "\$BC_TMP" ;;
  esac

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
