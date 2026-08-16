#!/usr/bin/env bash
# Generate the whole-genome synteny dot plot (hybrid assembly vs. reference)
# for a strain.
# Run from WSL: bash tools/plot_assembly_dotplot.sh <strain> <display_name>
# Example:      bash tools/plot_assembly_dotplot.sh AB30 AB030
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: plot_assembly_dotplot.sh <strain> <display_name> [ori_ter_json]}"
DISPLAY_NAME="${2:?Usage: plot_assembly_dotplot.sh <strain> <display_name> [ori_ter_json]}"
# Strain IDs in analysis/assembly_comparison/ (e.g. "Lac4") don't always match
# mapping/reference/ dirnames (e.g. "Lac-4") -- pass an override if they differ.
ORI_TER_JSON="${3:-mapping/reference/$STRAIN/ori_ter.json}"

# run_logged.py itself needs python-docx, which lives in seqqc, not rplots;
# call it via the seqqc interpreter directly while Rscript runs from rplots' PATH.
"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate whole-genome synteny dot plot for $STRAIN: hybrid assembly vs. reference genome, from NucDiff/MUMmer alignment blocks, colored by orientation" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/assembly_comparison/plot_assembly_dotplot.R $STRAIN \"$DISPLAY_NAME\" analysis/assembly_comparison/$STRAIN $ORI_TER_JSON" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/assembly_comparison/plot_assembly_dotplot.R \
     "$STRAIN" "$DISPLAY_NAME" "analysis/assembly_comparison/$STRAIN" "$ORI_TER_JSON"
