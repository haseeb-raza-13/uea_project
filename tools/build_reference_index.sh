#!/usr/bin/env bash
# Build a Bowtie2 index for a strain's reference genome.
# Run from WSL: bash tools/build_reference_index.sh <strain> <reference.fna>
# Example:      bash tools/build_reference_index.sh AB30 short_read_data/short_read_seq_1/AB30.fna
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

STRAIN="${1:?Usage: build_reference_index.sh <strain> <reference.fna>}"
REFERENCE="${2:?Usage: build_reference_index.sh <strain> <reference.fna>}"
OUT_DIR="mapping/reference/$STRAIN"
INDEX_PREFIX="$OUT_DIR/${STRAIN}_index"

mkdir -p "$OUT_DIR"
cp -n "$REFERENCE" "$OUT_DIR/$(basename "$REFERENCE")"

python3 scripts/run_logged.py \
  --purpose "Build Bowtie2 index for the $STRAIN reference genome" \
  --label "bowtie2-build $REFERENCE $INDEX_PREFIX" \
  -- bowtie2-build "$REFERENCE" "$INDEX_PREFIX"
