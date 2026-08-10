#!/usr/bin/env bash
# Align + deduplicate all 8 mff samples against the mff reference.
# Run from WSL: bash tools/run_mapping.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"
mkdir -p mapping/bam/mff

python3 scripts/run_logged.py \
  --purpose "Align mff short-read samples to the mff reference (Bowtie2) and remove PCR duplicates (samtools), required before marker-frequency analysis" \
  --tool bowtie2 \
  --label "bowtie2 + samtools sort/fixmate/markdup on 8 mff paired-end samples (PID-2861-25..32) -> mapping/bam/mff" \
  -- python3 tools/run_mapping.py
