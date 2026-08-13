#!/usr/bin/env bash
# Generate docs/FoldChange_Interpretation.docx (not committed; docs/ is untracked).
# Run from WSL: bash tools/write_fc_interpretation_report.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

python3 scripts/run_logged.py \
  --purpose "Write interpretation report for the Meropenem-vs-untreated log2 fold-change analysis (pooled test, transposed version, and per-window volcano plot) across all 5 strains" \
  --tool python3 \
  --label "write_fc_interpretation_report.py -> docs/FoldChange_Interpretation.docx (Times New Roman, headings 14pt, body 12pt)" \
  -- python3 tools/write_fc_interpretation_report.py
