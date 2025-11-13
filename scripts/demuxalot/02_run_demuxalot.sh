#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: 02_run_demuxalot.sh --sample ID --bam BAM --vcf VCF --barcodes FILE --inds FILE --outdir DIR [options]

Wrapper around Demuxalot.py inside the Demuxafy Singularity image.

Required:
  --sample ID         Pool/sample identifier (used in log messages)
  --bam PATH          Coordinate-sorted BAM/CRAM with barcode/UMI tags
  --vcf PATH          BAM-ordered VCF/BCF (bgzip + index, preferably filtered)
  --barcodes FILE     Barcode whitelist (one per line)
  --inds FILE         Donor list (matches sample IDs in VCF)
  --outdir DIR        Output directory for Demuxalot results

Options:
  --threads INT       Number of threads to pass via -p (default: 16)
  --cell-tag STR      Barcode SAM tag (default: CB)
  --umi-tag STR       UMI SAM tag (default: UB)
  --refine BOOL       Run genotype refinement? (default: true)
  --container PATH    Path to Demuxafy Singularity image (default: \$DEMUXAFY_SIF)
  --host-root PATH    Host directory to bind inside Singularity (default: /home/pr422)
  -h, --help          Show this message

Example:
  export DEMUXAFY_SIF=/path/to/Demuxafy.sif
  bash 02_run_demuxalot.sh --sample S2A --bam S2A.bam --vcf S2A.vcf.gz \\
       --barcodes S2A.barcodes.txt --inds S2A.inds.txt --outdir /tmp/S2A_demuxalot
USAGE
}

SAMPLE=""
BAM=""
VCF=""
BARCODES=""
INDS=""
OUTDIR=""
THREADS=16
CELL_TAG=CB
UMI_TAG=UB
REFINE=true
CONTAINER=${DEMUXAFY_SIF:-}
HOST_ROOT=${HOST_ROOT:-/home/pr422}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --bam) BAM="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --barcodes) BARCODES="$2"; shift 2 ;;
    --inds) INDS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --cell-tag) CELL_TAG="$2"; shift 2 ;;
    --umi-tag) UMI_TAG="$2"; shift 2 ;;
    --refine) REFINE="$2"; shift 2 ;;
    --container) CONTAINER="$2"; shift 2 ;;
    --host-root) HOST_ROOT="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "[ERROR] Unknown arg: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

for var in SAMPLE BAM VCF BARCODES INDS OUTDIR; do
  if [[ -z "${!var}" ]]; then
    echo "[ERROR] --${var,,} is required" >&2
    show_help >&2
    exit 1
  fi
done

for f in "$BAM" "$VCF" "$BARCODES" "$INDS"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Missing file: $f" >&2
    exit 1
  fi
done

if [[ ! -f "$BAM.bai" && ! -f "${BAM%.bam}.bai" && ! -f "$BAM.crai" ]]; then
  echo "[ERROR] BAM/CRAM index not found for $BAM" >&2
  exit 1
fi

if [[ ! -f "$VCF.tbi" && ! -f "$VCF.csi" ]]; then
  echo "[ERROR] VCF index (.tbi/.csi) is required for $VCF" >&2
  exit 1
fi

if [[ -z "$CONTAINER" || ! -f "$CONTAINER" ]]; then
  echo "[ERROR] Set --container or export DEMUXAFY_SIF to point to Demuxafy.sif" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

REFINE_ARG="True"
case "${REFINE,,}" in
  true|1|yes) REFINE_ARG=True ;;
  false|0|no) REFINE_ARG=False ;;
  *) echo "[ERROR] --refine must be true/false (got $REFINE)" >&2; exit 1 ;;
esac

CMD=(Demuxalot.py
  -b "$BARCODES"
  -a "$BAM"
  -n "$INDS"
  -v "$VCF"
  -o "$OUTDIR"
  -p "$THREADS"
  -r "$REFINE_ARG"
)

[[ -n "$CELL_TAG" ]] && CMD+=(-c "$CELL_TAG")
[[ -n "$UMI_TAG" ]] && CMD+=(-u "$UMI_TAG")

echo "[INFO] Running Demuxalot for $SAMPLE -> $OUTDIR"
singularity exec --bind "$HOST_ROOT":"$HOST_ROOT" "$CONTAINER" "${CMD[@]}"
echo "[INFO] Demuxalot completed for $SAMPLE"
