#!/usr/bin/env bash
# Generate the 6 publication-ready replication-profile plots for mff.
# Run from WSL: bash tools/plot_replication_profile.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

# run_logged.py itself needs python-docx, which lives in seqqc, not rplots;
# call it via the seqqc interpreter directly while Rscript runs from rplots' PATH.
"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Generate the 6 publication-ready replication-profile plots for mff (Skovgaard et al. 2011 methodology): A vs B and drug vs ND, for stationary and exponential phase" \
  --tool "$HOME/miniforge3/envs/rplots/bin/R" \
  --label "Rscript analysis/replication_profile/mff/plot_replication_profile.R mapping/coverage/mff mapping/coverage/sample_metadata.csv mff" \
  -- "$HOME/miniforge3/envs/rplots/bin/Rscript" analysis/replication_profile/mff/plot_replication_profile.R \
     mapping/coverage/mff mapping/coverage/sample_metadata.csv mff
