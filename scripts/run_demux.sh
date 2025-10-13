#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   run_demux.sh \
#     --bam /data/S10A/outs/possorted_genome_bam.bam \
#     --barcodes /data/S10A/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
#     --vcf-pileup /vcf/healthy_all.for_pileup.vcf.gz \
#     --vcf-demux  /vcf/healthy_all.merged.biallelic.vcf.gz \
#     --outdir /out \
#     --prefix S10A_mapped \
#     [--sm-list /lists/S10A_donors.txt] \
#     [--field GT|GP|GL] \
#     [--doublet-prior 0.1] \
#     [--tag-group CB] [--tag-UMI UB]
#
# Notes:
# - For 10x, default tags are CB (cell barcode) and UB (UMI). :contentReference[oaicite:4]{index=4}
# - The pileup VCF must carry AC and AN (use bcftools +fill-tags in Step 4.1). :contentReference[oaicite:5]{index=5}
# - Use --field GP for imputed genotypes; GT for hard-called genotypes. (popscle supports GT/GP/GL.) :contentReference[oaicite:6]{index=6}

# ---------- parse args ----------
BAM=""
BARCODES=""
VCF_PILEUP=""
VCF_DEMUX=""
OUTDIR=""
PREFIX=""
SM_LIST=""
FIELD="GT"
DOUBLETPRIOR="0.10"
TAG_GROUP="CB"
TAG_UMI="UB"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2;;
    --barcodes) BARCODES="$2"; shift 2;;
    --vcf-pileup) VCF_PILEUP="$2"; shift 2;;
    --vcf-demux) VCF_DEMUX="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --prefix) PREFIX="$2"; shift 2;;
    --sm-list) SM_LIST="$2"; shift 2;;
    --field) FIELD="$2"; shift 2;;
    --doublet-prior) DOUBLETPRIOR="$2"; shift 2;;
    --tag-group) TAG_GROUP="$2"; shift 2;;
    --tag-UMI) TAG_UMI="$2"; shift 2;;
    *) echo "[ERROR] Unknown arg: $1"; exit 1;;
  esac
done

for v in BAM BARCODES VCF_PILEUP VCF_DEMUX OUTDIR PREFIX; do
  if [[ -z "${!v}" ]]; then echo "[ERROR] Missing --${v//_/-}"; exit 1; fi
done

mkdir -p "$OUTDIR"
PLP_PREFIX="$OUTDIR/${PREFIX}_pileup"

# Decompress barcodes if needed (dsc-pileup expects a plain list)
BC_TMP="$BARCODES"
if [[ "$BARCODES" == *.gz ]]; then
  BC_TMP="$(mktemp)"
  gzip -cd "$BARCODES" > "$BC_TMP"
fi

echo "[INFO] Running dsc-pileup..."
popscle dsc-pileup \
  --sam "$BAM" \
  --vcf "$VCF_PILEUP" \
  --group-list "$BC_TMP" \
  --tag-group "$TAG_GROUP" \
  --tag-UMI "$TAG_UMI" \
  --out "$PLP_PREFIX"

echo "[INFO] Running demuxlet..."
DEMUX_DIR="$OUTDIR/${PREFIX}_demuxlet"
mkdir -p "$DEMUX_DIR"
set -x
popscle demuxlet \
  --plp "$PLP_PREFIX" \
  --vcf "$VCF_DEMUX" \
  --field "$FIELD" \
  ${SM_LIST:+--sm-list "$SM_LIST"} \
  --doublet-prior "$DOUBLETPRIOR" \
  --out "$DEMUX_DIR/demuxlet"
set +x

echo "[INFO] Done. Key outputs:"
echo "  - $DEMUX_DIR/demuxlet.best  (barcode → donor assignment, doublets, posteriors)"
echo "  - ${PLP_PREFIX}.{plp,var,cel}  (reusable pileups)"
