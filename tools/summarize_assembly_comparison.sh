#!/usr/bin/env bash
# Consolidate FastANI + NucDiff assembly-comparison output into one cross-strain
# summary table. Run this AFTER scripts/compare_assemblies.sh has completed on
# the HPC cluster and analysis/assembly_comparison/ has been copied back here.
# Run from WSL: bash tools/summarize_assembly_comparison.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

python3 scripts/run_logged.py \
  --purpose "Consolidate per-strain FastANI ANI/fragment-mapping stats (and, once parsed, NucDiff SNP/rearrangement counts) into one cross-strain assembly-comparison summary table" \
  --tool python3 \
  --label "summarize_assembly_comparison.py -> analysis/assembly_comparison/summary_table.csv" \
  -- python3 tools/summarize_assembly_comparison.py
