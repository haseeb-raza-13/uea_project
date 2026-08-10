#!/usr/bin/env bash
# Generate docs/mff_Replication_Profile_Findings.docx.
# Run from WSL: bash tools/write_mff_findings_report.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate seqqc

cd "$PROJECT_ROOT"

python3 scripts/run_logged.py \
  --purpose "Write formatted findings/discussion report for the mff replication-profile analysis (reasons for the flat signal, strategies, chromosomal-rearrangement and long-read-sequencing recommendations)" \
  --tool python3 \
  --label "write_mff_findings_report.py -> docs/mff_Replication_Profile_Findings.docx (Times New Roman, headings 14pt, body 12pt)" \
  -- python3 tools/write_mff_findings_report.py
