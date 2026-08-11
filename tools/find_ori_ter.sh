#!/usr/bin/env bash
# Independently locate oriC/ter on a strain's reference via cumulative GC-skew.
# Run from WSL: bash tools/find_ori_ter.sh <strain> <reference.fna> [approx_oric]
# Example:      bash tools/find_ori_ter.sh AB30 short_read_data/short_read_seq_1/AB30.fna
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: find_ori_ter.sh <strain> <reference.fna> [approx_oric]}"
REFERENCE="${2:?Usage: find_ori_ter.sh <strain> <reference.fna> [approx_oric]}"
APPROX_ORIC="${3:-}"
OUT_DIR="mapping/reference/$STRAIN"

mkdir -p "$OUT_DIR"

EXTRA_ARGS=()
LABEL_SUFFIX=""
if [ -n "$APPROX_ORIC" ]; then
  EXTRA_ARGS=(--approx-oric "$APPROX_ORIC")
  LABEL_SUFFIX=" --approx-oric $APPROX_ORIC"
fi

python3 scripts/run_logged.py \
  --purpose "Locate oriC/ter on the $STRAIN reference via cumulative GC-skew (Lobry's method)" \
  --tool python3 \
  --label "find_ori_ter.py $(basename "$REFERENCE") --window 1000$LABEL_SUFFIX" \
  -- python3 tools/find_ori_ter.py "$REFERENCE" "$OUT_DIR" --window 1000 "${EXTRA_ARGS[@]}"
