process HARMONIZE {
  tag "${sample_id}"
  label 'HARMONIZE'

  input:
  tuple val(sample_id), path(bam), path(vcf_in)

  output:
  tuple val(sample_id), path("${sample_id}.bamlike.af.clean.vcf.gz")

  script:
  """
  set -euo pipefail

  # 1) reorder header/contigs like BAM
  harmonize_vcf_bam.sh "${bam}" "${vcf_in}" z > ${sample_id}.bamlike.vcf.gz
  tabix -f -p vcf ${sample_id}.bamlike.vcf.gz

  # 2) add AF/AC/AN
  bcftools +fill-tags ${sample_id}.bamlike.vcf.gz -Oz -o ${sample_id}.bamlike.af.vcf.gz -- -t AC,AN,AF
  tabix -f -p vcf ${sample_id}.bamlike.af.vcf.gz

  # 3) keep only polymorphic biallelic SNPs (no no-ALT records)
  bcftools view -v snps -m2 -M2 -i 'N_ALT=1 && INFO/AC>0' \
    -Oz -o ${sample_id}.bamlike.af.clean.vcf.gz ${sample_id}.bamlike.af.vcf.gz
  tabix -f -p vcf ${sample_id}.bamlike.af.clean.vcf.gz
  """
}
