#!/usr/bin/env bash
# Remove Illumina/Nextera adapter contamination from all short-read samples (BBDuk).
# Run from WSL: bash tools/run_bbduk_trim.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"
mkdir -p trimmed_data/short_read/bbduk_stats

python3 scripts/run_logged.py \
  --purpose "Remove Illumina/Nextera adapter contamination from short-read data (BBDuk)" \
  --tool bbduk.sh \
  --label "bbduk.sh ktrim=r k=23 mink=11 hdist=1 tpe tbo ref=adapters.fa on 40 paired-end samples (short_read_seq_1/2/3 -> trimmed_data/short_read)" \
  -- python3 tools/run_bbduk_trim.py
