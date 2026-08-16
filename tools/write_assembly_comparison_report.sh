#!/usr/bin/env bash
# Generate docs/Assembly_Comparison_Report.docx (not committed; docs/ is untracked).
# Run from WSL: bash tools/write_assembly_comparison_report.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

"$HOME/miniforge3/envs/seqqc/bin/python3" tools/write_assembly_comparison_report.py
