#!/usr/bin/env bash
# Compute windowed, oriC/ter-normalized marker-frequency TSVs for a strain.
# Run from WSL: bash tools/compute_marker_frequency.sh <strain>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: compute_marker_frequency.sh <strain>}"
mkdir -p "mapping/coverage/$STRAIN"

python3 scripts/run_logged.py \
  --purpose "Compute 1kb/10kb windowed marker-frequency TSVs (read-start counts, log2 ratio to genome median, oriC/ter-normalized position) for all $STRAIN samples" \
  --tool bedtools \
  --label "bedtools genomecov -5 -bg + windowed binning (1kb, 10kb) on $STRAIN dedup BAMs -> mapping/coverage/$STRAIN" \
  -- python3 tools/compute_marker_frequency.py "$STRAIN"
