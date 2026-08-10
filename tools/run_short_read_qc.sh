#!/usr/bin/env bash
# Batch read QC (FastQC + MultiQC) for short-read FASTQ files. Run from WSL:
#   bash tools/run_short_read_qc.sh [--input-dir DIR ...] [--out-name NAME] [threads]
#
# With no arguments, QCs the raw data (short_read_data/short_read_seq_1/2/3) into
# qc_output/fastqc + qc_output/multiqc, exactly as before.
# --input-dir can be repeated to scan other directories (e.g. trimmed_data/short_read).
# --out-name NAME routes output to qc_output/fastqc_NAME + qc_output/multiqc_NAME.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

SHORT_READ_DIR="$PROJECT_ROOT/short_read_data"
SUBDIRS=()
OUT_NAME=""

while [ $# -gt 0 ]; do
  case "$1" in
    --input-dir)
      SUBDIRS+=("$2")
      shift 2
      ;;
    --out-name)
      OUT_NAME="$2"
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

if [ "${#SUBDIRS[@]}" -eq 0 ]; then
  SUBDIRS=("$SHORT_READ_DIR/short_read_seq_1" "$SHORT_READ_DIR/short_read_seq_2" "$SHORT_READ_DIR/short_read_seq_3")
fi

if [ -n "$OUT_NAME" ]; then
  FASTQC_OUT="$PROJECT_ROOT/qc_output/fastqc_$OUT_NAME"
  MULTIQC_OUT="$PROJECT_ROOT/qc_output/multiqc_$OUT_NAME"
  PURPOSE_SUFFIX=" ($OUT_NAME data)"
else
  FASTQC_OUT="$PROJECT_ROOT/qc_output/fastqc"
  MULTIQC_OUT="$PROJECT_ROOT/qc_output/multiqc"
  PURPOSE_SUFFIX=""
fi

mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

mapfile -t FASTQ_FILES < <(find "${SUBDIRS[@]}" -maxdepth 1 \( -iname "*.fastq" -o -iname "*.fastq.gz" \) | sort)
NUM_FILES=${#FASTQ_FILES[@]}

if [ "$NUM_FILES" -eq 0 ]; then
  echo "No FASTQ files found under: ${SUBDIRS[*]} — nothing to do." >&2
  exit 1
fi

if [ "$NUM_FILES" -ne 80 ]; then
  echo "Warning: expected 80 FASTQ files (40 samples x R1/R2), found $NUM_FILES." >&2
fi

THREADS="${1:-$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))}"
if [ "$THREADS" -gt "$NUM_FILES" ]; then
  THREADS=$NUM_FILES
fi

SUBDIR_BASENAMES=()
for d in "${SUBDIRS[@]}"; do SUBDIR_BASENAMES+=("$(basename "$d")"); done
SUBDIR_NAMES="$(printf '%s, ' "${SUBDIR_BASENAMES[@]}" | sed 's/, $//')"
echo "Found $NUM_FILES FASTQ files. Running FastQC with $THREADS thread(s)..."

python3 "$SCRIPT_DIR/../scripts/run_logged.py" \
  --purpose "Read quality control (FastQC) on all FASTQ files${PURPOSE_SUFFIX}" \
  --label "fastqc -t $THREADS -o ${FASTQC_OUT#$PROJECT_ROOT/} <$NUM_FILES FASTQ files from $SUBDIR_NAMES>" \
  -- fastqc -t "$THREADS" -o "$FASTQC_OUT" "${FASTQ_FILES[@]}"

echo "FastQC done. Running MultiQC to aggregate reports..."

python3 "$SCRIPT_DIR/../scripts/run_logged.py" \
  --purpose "Aggregate all FastQC reports into a single MultiQC summary${PURPOSE_SUFFIX}" \
  -- multiqc "$FASTQC_OUT" -o "$MULTIQC_OUT"

echo
echo "Done."
echo "FastQC reports: $FASTQC_OUT"
echo "MultiQC summary: $MULTIQC_OUT/multiqc_report.html"
