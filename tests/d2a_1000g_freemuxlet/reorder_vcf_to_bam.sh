#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 3 ]; then
  echo "Usage: $0 input.vcf.gz input.bam output.vcf.gz" >&2
  exit 1
fi

IN_VCF="$1"
IN_BAM="$2"
OUT_VCF="$3"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

> "$TMP_DIR/contigs.tsv"
while IFS=$'\t' read -r -a fields; do
  if [[ ${fields[0]} == "@SQ" ]]; then
    contig=""
    length=""
    for field in "${fields[@]:1}"; do
      case $field in
        SN:*) contig=${field#SN:} ;;
        LN:*) length=${field#LN:} ;;
      esac
    done
    if [[ -n $contig ]]; then
      [[ -z $length ]] && length=1
      printf '%s\t%s\n' "$contig" "$length" >> "$TMP_DIR/contigs.tsv"
    fi
  fi
done < <(samtools view -H "$IN_BAM")

awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' "$TMP_DIR/contigs.tsv" > "$TMP_DIR/contigs_header.txt"

bcftools view -h "$IN_VCF" > "$TMP_DIR/orig_header.txt"

awk 'NR==FNR {contig_lines[++n] = $0; next} /^##contig=/ {next} /^#CHROM/ {for (i = 1; i <= n; i++) print contig_lines[i]; print $0; next} {print}' \
  "$TMP_DIR/contigs_header.txt" "$TMP_DIR/orig_header.txt" > "$TMP_DIR/new_header.txt"

bcftools index -f "$IN_VCF"

cp "$IN_VCF" "$TMP_DIR/input.vcf.gz"
cp "$IN_VCF".csi "$TMP_DIR/input.vcf.gz.csi" 2>/dev/null || true
cp "$IN_VCF".tbi "$TMP_DIR/input.vcf.gz.tbi" 2>/dev/null || true

bcftools reheader -h "$TMP_DIR/new_header.txt" "$TMP_DIR/input.vcf.gz" -o "$TMP_DIR/reheader.vcf.gz"

tabix -f -p vcf "$TMP_DIR/reheader.vcf.gz"

bcftools sort -O z -o "$OUT_VCF" "$TMP_DIR/reheader.vcf.gz"

tabix -f -p vcf "$OUT_VCF"
