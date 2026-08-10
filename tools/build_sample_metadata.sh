#!/usr/bin/env bash
# Parse Sample Name illumina.xlsx into mapping/coverage/sample_metadata.csv.
# Run from WSL: bash tools/build_sample_metadata.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"
mkdir -p mapping/coverage

python3 scripts/run_logged.py \
  --purpose "Parse Sample Name illumina.xlsx into a tidy per-sample metadata CSV (strain, replicate, treatment, phase) for joining against marker-frequency TSVs" \
  -- python3 tools/build_sample_metadata.py
