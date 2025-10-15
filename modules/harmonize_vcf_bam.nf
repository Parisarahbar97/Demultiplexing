process harmonize_vcf_bam {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/harmonize", mode: 'copy'

  input:
  tuple val(sample_id), path(bam), path(vcf)

  output:
  tuple val(sample_id),
        path("${sample_id}.bamlike.vcf.gz"),
        path("${sample_id}.bamlike.vcf.gz.tbi"),
        path("${sample_id}.bamlike.af.vcf.gz"),
        path("${sample_id}.bamlike.af.vcf.gz.tbi"),
        path("${sample_id}.bamlike.af.clean.vcf.gz"),
        path("${sample_id}.bamlike.af.clean.vcf.gz.tbi")

  script:
  def filterExpr = params.vcf_filter == 'biallelic_all'
      ? '-i "N_ALT=1"'
      : '-i "N_ALT=1 && INFO/AC>0"'
  """
  set -euo pipefail

  # 1) Harmonize contig order to match BAM
  harmonize_vcf_bam.sh "$bam" "$vcf" z > ${sample_id}.bamlike.vcf.gz
  tabix -f -p vcf ${sample_id}.bamlike.vcf.gz

  # 2) Fill AC/AN/AF
  bcftools +fill-tags ${sample_id}.bamlike.vcf.gz -Oz -o ${sample_id}.bamlike.af.vcf.gz -- -t AC,AN,AF
  tabix -f -p vcf ${sample_id}.bamlike.af.vcf.gz

  # 3) Clean to biallelic SNPs; optionally restrict to polymorphic
  bcftools view -v snps -m2 -M2 ${filterExpr} \
    -Oz -o ${sample_id}.bamlike.af.clean.vcf.gz ${sample_id}.bamlike.af.vcf.gz
  tabix -f -p vcf ${sample_id}.bamlike.af.clean.vcf.gz
  """
}
