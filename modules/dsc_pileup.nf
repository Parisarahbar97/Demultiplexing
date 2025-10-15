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
  """
  set -euo pipefail
  BC_TMP=\$(mktemp)
  case "$barcodes" in
    *.gz) gzip -cd "$barcodes" > "\$BC_TMP" ;;
    *)    cp "$barcodes" "\$BC_TMP" ;;
  esac

  popscle dsc-pileup \
    --sam "$bam" \
    --vcf "$vcf_for_pileup" \
    --group-list "\$BC_TMP" \
    --tag-group CB --tag-UMI UB \
    --out "${sample_id}_pileup"

  rm -f "\$BC_TMP"
  """
}
