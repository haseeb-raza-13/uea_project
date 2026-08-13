#!/usr/bin/env bash
# Generate the fold-change/volcano-plot replication-profile figures for a strain:
# Panels A/B are identical to the log2 version; Panel C is a classic volcano
# plot, log2(Meropenem/untreated) per 1kb window (x) vs -log10(per-window
# t-test p-value) (y), with Mann-Whitney U p-values also computed and written
# to the results CSV.
# Run from WSL: bash tools/plot_replication_profile_fc_volcano.sh <strain> <display_name>
# Example:      bash tools/plot_replication_profile_fc_volcano.sh AB30 AB030
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: plot_replication_profile_fc_volcano.sh <strain> <display_name>}"
DISPLAY_NAME="${2:?Usage: plot_replication_profile_fc_volcano.sh <strain> <display_name>}"

"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate fold-change/volcano-plot replication-profile figures for $STRAIN: Panel C = log2(Meropenem/untreated) per 1kb window vs -log10(per-window t-test p), Mann-Whitney U p-values also computed, for stationary and exponential phase" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/replication_profile/plot_replication_profile_fc_volcano.R $STRAIN \"$DISPLAY_NAME\" mapping/coverage/$STRAIN mapping/coverage/sample_metadata.csv analysis/replication_profile/$STRAIN/figures_fc_volcano" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/replication_profile/plot_replication_profile_fc_volcano.R \
     "$STRAIN" "$DISPLAY_NAME" "mapping/coverage/$STRAIN" "mapping/coverage/sample_metadata.csv" "analysis/replication_profile/$STRAIN/figures_fc_volcano"
