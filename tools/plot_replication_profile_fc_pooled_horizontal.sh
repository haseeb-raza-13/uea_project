#!/usr/bin/env bash
# Generate the transposed pooled-test fold-change replication-profile figures
# for a strain: Panels A/B are identical to the log2 version; Panel C shows
# the same log2(Meropenem/untreated) per 1kb window and the same overall
# t-test / Mann-Whitney U p-values as figures_fc_pooled, but with axes rotated
# -- genomic position on x, log2FC on y, with solid horizontal +-1 threshold
# lines instead of dashed vertical ones.
# Run from WSL: bash tools/plot_replication_profile_fc_pooled_horizontal.sh <strain> <display_name>
# Example:      bash tools/plot_replication_profile_fc_pooled_horizontal.sh AB30 AB030
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: plot_replication_profile_fc_pooled_horizontal.sh <strain> <display_name>}"
DISPLAY_NAME="${2:?Usage: plot_replication_profile_fc_pooled_horizontal.sh <strain> <display_name>}"

"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate transposed pooled-test fold-change replication-profile figures for $STRAIN: Panel C = log2(Meropenem/untreated) per 1kb window (y) vs genomic position (x), solid horizontal +-1 threshold lines, same overall t-test/Mann-Whitney U p-values as figures_fc_pooled, for stationary and exponential phase" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/replication_profile/plot_replication_profile_fc_pooled_horizontal.R $STRAIN \"$DISPLAY_NAME\" mapping/coverage/$STRAIN mapping/coverage/sample_metadata.csv analysis/replication_profile/$STRAIN/figures_fc_pooled_horizontal" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/replication_profile/plot_replication_profile_fc_pooled_horizontal.R \
     "$STRAIN" "$DISPLAY_NAME" "mapping/coverage/$STRAIN" "mapping/coverage/sample_metadata.csv" "analysis/replication_profile/$STRAIN/figures_fc_pooled_horizontal"
