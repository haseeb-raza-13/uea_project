#!/usr/bin/env bash
# Generate the 2 publication-ready composite replication-profile figures
# (stationary phase, exponential phase) for a strain.
# Run from WSL: bash tools/plot_replication_profile.sh <strain> <display_name>
# Example:      bash tools/plot_replication_profile.sh AB30 AB030
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: plot_replication_profile.sh <strain> <display_name>}"
DISPLAY_NAME="${2:?Usage: plot_replication_profile.sh <strain> <display_name>}"

# run_logged.py itself needs python-docx, which lives in seqqc, not rplots;
# call it via the seqqc interpreter directly while Rscript runs from rplots' PATH.
"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate publication-ready composite replication-profile figures for $STRAIN (Skovgaard et al. 2011 methodology): drug-treated and untreated, replicate A vs B and treatment comparison, for stationary and exponential phase" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/replication_profile/plot_replication_profile.R $STRAIN \"$DISPLAY_NAME\" mapping/coverage/$STRAIN mapping/coverage/sample_metadata.csv analysis/replication_profile/$STRAIN/figures" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/replication_profile/plot_replication_profile.R \
     "$STRAIN" "$DISPLAY_NAME" "mapping/coverage/$STRAIN" "mapping/coverage/sample_metadata.csv" "analysis/replication_profile/$STRAIN/figures"
