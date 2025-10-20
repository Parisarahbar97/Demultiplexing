#!/usr/bin/env bash
#
# run_plink_qc.sh
# ----------------
# Per-pool pre-imputation quality control that mirrors the manual pipeline
# used during development. This script expects PLINK 2, PLINK 1.9 and
# bcftools (>=1.10) on PATH – the `Docker/Dockerfile.impute` image provides
# these tools.
#
# Usage:
#   ./scripts/run_plink_qc.sh \
#       --manifest imputation_work/00_raw/vcf_manifest.txt \
#       --outdir   imputation_work/01_plink_qc \
#       --logdir   imputation_work/logs \
#       --threads  8
#
# The manifest must contain one absolute VCF path per line. The script will
# create PLINK PGENs, basic QC metrics and cleaned VCFs suitable for the
# TOPMed preparation step.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_plink_qc.sh [options]

Required options:
  -m, --manifest FILE   Text file with one VCF (.vcf or .vcf.gz) per line.
  -o, --outdir   DIR    Directory for PLINK/PVCF outputs.
  -l, --logdir   DIR    Directory for PLINK/BCFtools logs.

Optional:
  -t, --threads INT     Threads to pass to PLINK2 (default: 8).
  -h, --help            Show this help and exit.
USAGE
}

manifest=""
outdir=""
logdir=""
threads=8

while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--manifest) manifest="$2"; shift 2 ;;
    -o|--outdir)   outdir="$2";   shift 2 ;;
    -l|--logdir)   logdir="$2";   shift 2 ;;
    -t|--threads)  threads="$2";  shift 2 ;;
    -h|--help)     usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "${manifest}" ]] && { echo "ERROR: --manifest is required" >&2; usage; exit 1; }
[[ -z "${outdir}"   ]] && { echo "ERROR: --outdir is required" >&2; usage; exit 1; }
[[ -z "${logdir}"   ]] && { echo "ERROR: --logdir is required" >&2; usage; exit 1; }

manifest=$(realpath "${manifest}")
outdir=$(realpath -m "${outdir}")
logdir=$(realpath -m "${logdir}")

[[ -f "${manifest}" ]] || { echo "ERROR: manifest not found: ${manifest}" >&2; exit 1; }

mkdir -p "${outdir}" "${logdir}"

command -v plink2 >/dev/null 2>&1 || { echo "ERROR: plink2 not found on PATH" >&2; exit 1; }
command -v plink  >/dev/null 2>&1 || { echo "ERROR: plink (1.9) not found on PATH" >&2; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found on PATH" >&2; exit 1; }

mapfile -t vcfs < <(grep -v '^\s*$' "${manifest}")
if [[ "${#vcfs[@]}" -eq 0 ]]; then
  echo "ERROR: manifest is empty: ${manifest}" >&2
  exit 1
fi

echo "[run_plink_qc] Using manifest: ${manifest}"
echo "[run_plink_qc] Output directory: ${outdir}"
echo "[run_plink_qc] Logs directory:   ${logdir}"
echo "[run_plink_qc] Processing ${#vcfs[@]} VCF(s)…"

mapfile -t rename_chr < <(seq 1 22 | awk '{printf("chr%s\t%s\n",$1,$1)}')

for vcf in "${vcfs[@]}"; do
  vcf=$(realpath "${vcf}")
  if [[ ! -f "${vcf}" ]]; then
    echo "[run_plink_qc] WARNING: skipping missing file ${vcf}" >&2
    continue
  fi

  base=$(basename "${vcf}")
  prefix="${base%.vcf.gz}"
  prefix="${prefix%.vcf}"
  out_prefix="${outdir}/${prefix}"

  echo "[run_plink_qc] >>> ${prefix}"

  tmp_dir=$(mktemp -d)
  trap 'rm -rf "${tmp_dir}"' EXIT

  vcf_in="${vcf}"
  # If contigs already contain "chr", we do not need to rename.
  if bcftools view -h "${vcf}" | grep -q '^##contig=<ID=chr'; then
    :
  else
    mapfile -t tmp_map < <(printf "%s\n" "${rename_chr[@]}")
    chr_map="${tmp_dir}/chr_map.txt"
    printf "%s\n" "${tmp_map[@]}" > "${chr_map}"
    bcftools annotate --rename-chrs "${chr_map}" \
      -Oz -o "${tmp_dir}/${prefix}.chr.vcf.gz" "${vcf}"
    tabix -f -p vcf "${tmp_dir}/${prefix}.chr.vcf.gz"
    vcf_in="${tmp_dir}/${prefix}.chr.vcf.gz"
  fi

  # Main PLINK2 QC generating pgen/pvar/psam
  plink2 \
    --threads "${threads}" \
    --vcf "${vcf_in}" \
    --vcf-half-call missing \
    --set-missing-var-ids @:#:\$r:\$a \
    --chr-set 38 no-xy \
    --chr 1-22 \
    --snps-only just-acgt \
    --max-alleles 2 \
    --rm-dup force-first \
    --geno 0.05 \
    --hwe 1e-6 midp keep-fewhet \
    --make-pgen \
    --out "${out_prefix}" \
    >"${logdir}/${prefix}.plink2.qc.log" 2>&1

  # Metrics (allele counts + sample missingness)
  plink2 --threads "${threads}" \
    --pfile "${out_prefix}" \
    --freq counts \
    --missing sample-only \
    --out "${out_prefix}.qcmetrics" \
    >"${logdir}/${prefix}.metrics.log" 2>&1

  # Export cleaned, polymorphic, autosomal SNP VCF (MAC >=1)
  plink2 --threads "${threads}" \
    --pfile "${out_prefix}" \
    --chr 1-22 \
    --snps-only just-acgt \
    --max-alleles 2 \
    --mac 1 \
    --export vcf bgz id-paste=iid \
    --out "${out_prefix}.clean" \
    >"${logdir}/${prefix}.export.log" 2>&1

  # Drop PLINK's non-standard header line and ensure sorted output
  bcftools view -Ov "${out_prefix}.clean.vcf.gz" \
    | sed '/^##chrSet=/d' \
    | bcftools sort -Oz -o "${out_prefix}.clean.sorted.vcf.gz"

  tabix -f -p vcf "${out_prefix}.clean.sorted.vcf.gz"

  mv -f "${out_prefix}.clean.sorted.vcf.gz"     "${out_prefix}.clean.vcf.gz"
  mv -f "${out_prefix}.clean.sorted.vcf.gz.tbi" "${out_prefix}.clean.vcf.gz.tbi"

  # Clean tmp dir for next iteration
  rm -rf "${tmp_dir}"
  trap - EXIT
done

echo "[run_plink_qc] Completed. Outputs available under ${outdir}"
