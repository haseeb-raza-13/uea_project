#!/usr/bin/env bash
# Generate the circular genome-wide variant ideogram (SNP/indel density +
# structural rearrangements) for a strain.
# Run from WSL: bash tools/plot_assembly_circos.sh <strain> <display_name> [ori_ter_json]
# Example:      bash tools/plot_assembly_circos.sh AB30 AB030
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: plot_assembly_circos.sh <strain> <display_name> [ori_ter_json]}"
DISPLAY_NAME="${2:?Usage: plot_assembly_circos.sh <strain> <display_name> [ori_ter_json]}"
ORI_TER_JSON="${3:-mapping/reference/$STRAIN/ori_ter.json}"

"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate circular genome-wide variant ideogram for $STRAIN: SNP density, indel density, and structural rearrangement tracks from NucDiff output, oriC/ter marked" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/assembly_comparison/plot_assembly_circos.R $STRAIN \"$DISPLAY_NAME\" analysis/assembly_comparison/$STRAIN $ORI_TER_JSON" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/assembly_comparison/plot_assembly_circos.R \
     "$STRAIN" "$DISPLAY_NAME" "analysis/assembly_comparison/$STRAIN" "$ORI_TER_JSON"
