#!/usr/bin/env bash
set -euo pipefail

# Harmonize a VCF header/record order to match a BAM's contig order (for popscle).
# Optionally add INFO tags AC, AN, AF at the end.
#
# Usage:
#   harmonize_vcf_bam.sh BAM VCF_IN [--out-type v|z|u|b] [--add-af] [--af-tags 'AC,AN,AF'] [--tmpdir DIR]
#
# Notes:
#   - Writes the final VCF/BCF to STDOUT. Redirect it to a file.
#   - --out-type:
#       v = uncompressed VCF     z = bgzipped VCF
#       u = uncompressed BCF     b = compressed BCF
#     (default: z)
#   - --add-af will run `bcftools +fill-tags` with the tags in --af-tags (default AC,AN,AF)
#   - Requires: samtools, bcftools (with plugins), awk (mawk or gawk)

usage() {
  cat >&2 <<'EOF'
Usage:
  harmonize_vcf_bam.sh BAM VCF_IN [--out-type v|z|u|b] [--add-af] [--af-tags 'AC,AN,AF'] [--tmpdir DIR]

Examples:
  # Write bgzipped harmonized VCF to file:
  harmonize_vcf_bam.sh sample.bam pool.vcf z > pool.bamlike.vcf.gz

  # Same, and add AC/AN/AF in one go:
  harmonize_vcf_bam.sh sample.bam pool.vcf --out-type z --add-af > pool.bamlike.af.vcf.gz
EOF
}

# ---- parse args ----
if [[ $# -lt 2 ]]; then usage; exit 1; fi

BAM=""; VCF_IN=""
OUTTYPE="z"
ADD_AF=0
AF_TAGS="AC,AN,AF"
USER_TMPDIR=""

# positional BAM, VCF first
BAM="$1"; shift
VCF_IN="$1"; shift

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-type) OUTTYPE="${2:-z}"; shift 2;;
    --add-af)   ADD_AF=1; shift;;
    --af-tags)  AF_TAGS="${2:-AC,AN,AF}"; shift 2;;
    --tmpdir)   USER_TMPDIR="${2:-}"; shift 2;;
    -h|--help)  usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

# ---- setup temp ----
if [[ -n "$USER_TMPDIR" ]]; then
  tmpdir="$USER_TMPDIR"
  mkdir -p "$tmpdir"
else
  tmpdir="$(mktemp -d)"
  trap 'rm -rf "$tmpdir"' EXIT
fi

hdr_all="$tmpdir/hdr.all.txt"
meta_nocontig="$tmpdir/meta.nocontig.txt"
hdr_cols="$tmpdir/cols.txt"
vcf_contigs="$tmpdir/vcf.contigs.tsv"
bam_order="$tmpdir/bam.order.txt"
hdr_contigs_new="$tmpdir/hdr.contigs.new.txt"
hdr_new="$tmpdir/hdr.new.txt"

# ---- 1) Split VCF header into parts ----
bcftools view -h "$VCF_IN" > "$hdr_all"

# keep all header lines except contigs and #CHROM
grep -v '^##contig=<' "$hdr_all" | grep -v '^#CHROM' > "$meta_nocontig"
# the #CHROM line
grep '^#CHROM' "$hdr_all" > "$hdr_cols"

# map VCF contig -> length (safe for mawk)
awk -F'[=,>]' '
  /^##contig=<ID=/ {
    id=""; len="";
    for (i=1;i<=NF;i++) {
      if ($i=="ID")     { if (i<NF) id=$(i+1); }
      else if ($i=="length") { if (i<NF) len=$(i+1); }
    }
    if (id!="") { print id "\t" len; }
  }
' "$hdr_all" > "$vcf_contigs"

# ---- 2) Extract BAM contig order ----
samtools view -H "$BAM" \
| awk -F'\t' '
  $1=="@SQ" {
    id="";
    for(i=1;i<=NF;i++){
      if ($i ~ /^SN:/) { id=substr($i,4); }
    }
    if (id!="") { print id; }
  }
' > "$bam_order"

# ---- 3) Build new contig header in BAM order for contigs present in the VCF ----
awk '
  NR==FNR { len[$1]=$2; next }               # first file: contig -> length
  ($1 in len) {
    printf("##contig=<ID=%s,length=%s>\n",$1,len[$1]);
  }
' "$vcf_contigs" "$bam_order" > "$hdr_contigs_new"

# ---- 4) Assemble new header ----
cat "$meta_nocontig" "$hdr_contigs_new" "$hdr_cols" > "$hdr_new"

# ---- 5) Reheader + sort. Optionally add AF/AC/AN. Output to STDOUT. ----
if [[ "$ADD_AF" -eq 1 ]]; then
  bcftools reheader -h "$hdr_new" "$VCF_IN" \
  | bcftools sort -T "$tmpdir" -O v \
  | bcftools +fill-tags -O "$OUTTYPE" -- -t "$AF_TAGS"
else
  bcftools reheader -h "$hdr_new" "$VCF_IN" \
  | bcftools sort -T "$tmpdir" -O "$OUTTYPE"
fi
