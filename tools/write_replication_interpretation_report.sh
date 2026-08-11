#!/usr/bin/env bash
# Generate docs/Replication_Profile_Interpretation.docx.
# Run from WSL: bash tools/write_replication_interpretation_report.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

python3 scripts/run_logged.py \
  --purpose "Write formatted per-strain interpretation report for the replication-profile figures and Ori:Ter ratios across all 5 strains" \
  --tool python3 \
  --label "write_replication_interpretation_report.py -> docs/Replication_Profile_Interpretation.docx (Times New Roman, headings 14pt, body 12pt)" \
  -- python3 tools/write_replication_interpretation_report.py
