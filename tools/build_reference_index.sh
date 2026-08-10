#!/usr/bin/env bash
# Build a Bowtie2 index for the mff reference genome.
# Run from WSL: bash tools/build_reference_index.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

REFERENCE="short_read_data/short_read_seq_1/mff.fna"
OUT_DIR="mapping/reference/mff"
INDEX_PREFIX="$OUT_DIR/mff_index"

mkdir -p "$OUT_DIR"
cp -n "$REFERENCE" "$OUT_DIR/mff.fna"

python3 scripts/run_logged.py \
  --purpose "Build Bowtie2 index for the mff reference genome (CP012004.1)" \
  --label "bowtie2-build $REFERENCE $INDEX_PREFIX" \
  -- bowtie2-build "$REFERENCE" "$INDEX_PREFIX"
