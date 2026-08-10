#!/usr/bin/env bash
# Independently locate oriC/ter on the mff reference via cumulative GC-skew.
# Run from WSL: bash tools/find_ori_ter.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

REFERENCE="short_read_data/short_read_seq_1/mff.fna"
OUT_DIR="mapping/reference/mff"
APPROX_ORIC=1918991

mkdir -p "$OUT_DIR"

python3 scripts/run_logged.py \
  --purpose "Independently locate oriC/ter on mff.fna via cumulative GC-skew (Lobry's method), since the user-provided coordinates were tied to a chromosome length that didn't match this reference file" \
  --tool python3 \
  --label "find_ori_ter.py mff.fna --window 1000 --approx-oric $APPROX_ORIC" \
  -- python3 tools/find_ori_ter.py "$REFERENCE" "$OUT_DIR" --window 1000 --approx-oric "$APPROX_ORIC"
