#!/usr/bin/env bash
# Generate the cross-strain assembly-comparison summary composite (variant
# counts, ANI, assembly completeness) from analysis/assembly_comparison/summary_table.csv.
# Run from WSL: bash tools/plot_assembly_summary.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate cross-strain assembly-comparison summary composite: NucDiff variant-type counts, FastANI, and assembly completeness, for all 5 strains" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/assembly_comparison/plot_assembly_summary.R" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/assembly_comparison/plot_assembly_summary.R
