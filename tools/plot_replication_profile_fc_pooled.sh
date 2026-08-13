#!/usr/bin/env bash
# Generate the fold-change/pooled-test replication-profile figures for a strain:
# Panels A/B are identical to the log2 version; Panel C shows log2(Meropenem/
# untreated) per 1kb window vs genomic position, with ONE overall t-test and
# ONE overall Mann-Whitney U p-value pooling all windows as the sample.
# Run from WSL: bash tools/plot_replication_profile_fc_pooled.sh <strain> <display_name>
# Example:      bash tools/plot_replication_profile_fc_pooled.sh AB30 AB030
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: plot_replication_profile_fc_pooled.sh <strain> <display_name>}"
DISPLAY_NAME="${2:?Usage: plot_replication_profile_fc_pooled.sh <strain> <display_name>}"

"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate fold-change/pooled-test replication-profile figures for $STRAIN: Panel C = log2(Meropenem/untreated) per 1kb window vs genomic position, with an overall t-test and Mann-Whitney U p-value pooling all windows, for stationary and exponential phase" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/replication_profile/plot_replication_profile_fc_pooled.R $STRAIN \"$DISPLAY_NAME\" mapping/coverage/$STRAIN mapping/coverage/sample_metadata.csv analysis/replication_profile/$STRAIN/figures_fc_pooled" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/replication_profile/plot_replication_profile_fc_pooled.R \
     "$STRAIN" "$DISPLAY_NAME" "mapping/coverage/$STRAIN" "mapping/coverage/sample_metadata.csv" "analysis/replication_profile/$STRAIN/figures_fc_pooled"
