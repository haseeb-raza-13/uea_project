#!/usr/bin/env bash
# Generate the zoomed inversion-breakpoint diagrams for a strain (default:
# AB30, the only strain with NucDiff-called inversions).
# Run from WSL: bash tools/plot_assembly_locus.sh [strain] [display_name]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

STRAIN="${1:-AB30}"
DISPLAY_NAME="${2:-AB030}"

"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate zoomed inversion-breakpoint ribbon diagrams for $STRAIN: candidate inversion loci from NucDiff, with assembly-confidence caveats (contig identity, fragmentation)" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/assembly_comparison/plot_assembly_locus.R $STRAIN \"$DISPLAY_NAME\"" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/assembly_comparison/plot_assembly_locus.R \
     "$STRAIN" "$DISPLAY_NAME"
