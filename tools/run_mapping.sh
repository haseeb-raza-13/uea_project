#!/usr/bin/env bash
# Align + deduplicate all samples for a strain against its reference.
# Run from WSL: bash tools/run_mapping.sh <strain>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: run_mapping.sh <strain>}"
mkdir -p "mapping/bam/$STRAIN"

python3 scripts/run_logged.py \
  --purpose "Align $STRAIN short-read samples to the $STRAIN reference (Bowtie2) and remove PCR duplicates (samtools), required before marker-frequency analysis" \
  --tool bowtie2 \
  --label "bowtie2 + samtools sort/fixmate/markdup on $STRAIN paired-end samples -> mapping/bam/$STRAIN" \
  -- python3 tools/run_mapping.py "$STRAIN"
