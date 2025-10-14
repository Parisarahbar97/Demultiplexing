#!/usr/bin/env bash
set -euo pipefail
# Usage:
#   harmonize_vcf_bam.sh BAM.vbam VCF_IN.vcf[.gz] VCF_OUT.vcf.gz [--chr-map map.txt]
# 1) Get contig order from BAM -> rebuild VCF header in that order
# 2) Sort records to that order
# 3) Ensure INFO/AC, INFO/AN, INFO/AF exist (needed by dsc-pileup)

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 BAM VCF_IN VCF_OUT [--chr-map map.txt]" >&2; exit 1
fi
BAM="$1"; VCF_IN="$2"; VCF_OUT="$3"; CHR_MAP=""
shift 3 || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chr-map) CHR_MAP="$2"; shift 2;;
    *) echo "Unknown option: $1" >&2; exit 1;;
  esac
done

command -v bcftools >/dev/null || { echo "bcftools not found" >&2; exit 1; }
command -v samtools >/dev/null || { echo "samtools not found" >&2; exit 1; }
command -v tabix    >/dev/null || { echo "tabix not found" >&2; exit 1; }

workdir="$(mktemp -d)"; trap 'rm -rf "$workdir"' EXIT

# contigs from BAM in BAM order
samtools view -H "$BAM" \
| awk -F'\t' '$1=="@SQ"{sn="";ln="";
  for(i=2;i<=NF;i++){if($i~^"SN:")sn=substr($i,4);else if($i~^"LN:")ln=substr($i,4)}
  if(sn!="")printf "##contig=<ID=%s,length=%s>\n",sn,ln}' \
> "$workdir/contigs.from_bam.txt"

# optional rename of VCF contigs
VCF_STEP="$VCF_IN"
if [[ -n "$CHR_MAP" ]]; then
  VCF_STEP="$workdir/renamed.vcf.gz"
  bcftools annotate --rename-chrs "$CHR_MAP" "$VCF_IN" -Oz -o "$VCF_STEP"
  tabix -f -p vcf "$VCF_STEP"
fi

# rebuild header: meta (no contigs) + BAM contigs + #CHROM
bcftools view -h "$VCF_STEP" > "$workdir/hdr.all.txt"
grep -v '^##contig=<' "$workdir/hdr.all.txt" | grep -v '^#CHROM' > "$workdir/meta.nocontig.txt"
grep  '^#CHROM'       "$workdir/hdr.all.txt" > "$workdir/cols.txt"
cat "$workdir/meta.nocontig.txt" "$workdir/contigs.from_bam.txt" "$workdir/cols.txt" > "$workdir/hdr.bamorder.txt"

# apply header + sort records to that order
bamlike="$workdir/bamlike.vcf.gz"
bcftools reheader -h "$workdir/hdr.bamorder.txt" "$VCF_STEP" | bcftools sort -Oz -o "$bamlike"
tabix -f -p vcf "$bamlike"

# ensure AC/AN/AF present
bcftools +fill-tags "$bamlike" -Oz -o "$VCF_OUT" -- -t AC,AN,AF
tabix -f -p vcf "$VCF_OUT"

echo "Harmonized VCF: $VCF_OUT"
bcftools view -h "$VCF_OUT" | grep -E 'ID=(AC|AN|AF)' || true
