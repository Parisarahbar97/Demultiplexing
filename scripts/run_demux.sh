#!/usr/bin/env bash
set -euo pipefail

# Usage example:
#   run_demux.sh \
#     --bam /data/D1A_mapped/outs/possorted_genome_bam.bam \
#     --barcodes /data/D1A_mapped/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
#     --vcf-pileup /vcf/healthy_all.for_pileup.vcf.gz \
#     --vcf-demux  /vcf/healthy_all.merged.biallelic.vcf.gz \
#     --outdir /out \
#     --prefix D1A_mapped \
#     --sm-list /lists/D1A_pool_donors.txt \
#     --field GT

# Required:
BAM=""; BARCODES=""; VCF_PILEUP=""; VCF_DEMUX=""; OUTDIR=""; PREFIX=""
# Optional:
SM_LIST=""
FIELD="GT"          # GT for hard-called; GP for imputed (diseased)
DOUBLETPRIOR="0.10" # typical 5–10% expected doublets
TAG_GROUP="CB"      # 10x default
TAG_UMI="UB"        # 10x default

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2;;
    --barcodes) BARCODES="$2"; shift 2;;
    --vcf-pileup) VCF_PILEUP="$2"; shift 2;;
    --vcf-demux)  VCF_DEMUX="$2"; shift 2;;
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

# dsc-pileup expects a plain list; gunzip if needed
BC_TMP="$BARCODES"
if [[ "$BARCODES" == *.gz ]]; then
  BC_TMP="$(mktemp)"
  gzip -cd "$BARCODES" > "$BC_TMP"
fi

echo "[INFO] Running dsc-pileup..."
# popscle docs: 10x tags are CB (barcodes) and UB (UMIs); pileup VCF must have AC/AN. :contentReference[oaicite:1]{index=1}
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

echo "[INFO] Done:"
echo "  - $DEMUX_DIR/demuxlet.best   (barcode → donor; posteriors; doublet calls)"
echo "  - ${PLP_PREFIX}.{plp,var,cel}  (reusable pileups)"
