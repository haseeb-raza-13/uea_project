#!/usr/bin/env bash
# Compare each strain's hybrid assembly against its reference genome:
# FastANI (Average Nucleotide Identity) + NucDiff (SNPs / structural
# rearrangements). Runs entirely locally (no HPC/SLURM needed).
# Run from WSL: bash tools/compare_assemblies.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

source "$HOME/miniforge3/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "ani_nucdiff"; then
  echo "Environment 'ani_nucdiff' not found - creating with fastani and nucdiff (via mamba)."
  mamba create -y -n ani_nucdiff -c bioconda -c conda-forge fastani nucdiff
fi

cd "$PROJECT_ROOT"

# run_logged.py itself needs python-docx, which lives in seqqc, not
# ani_nucdiff; call it via the seqqc interpreter directly while
# compare_assemblies.py resolves fastANI/nucdiff from ani_nucdiff's bin/ dir.
"$HOME/miniforge3/envs/seqqc/bin/python3" scripts/run_logged.py \
  --purpose "Compare each strain's Unicycler hybrid assembly against its public reference genome: Average Nucleotide Identity (FastANI) and SNPs/structural rearrangements (NucDiff), for all 5 strains" \
  --tool "$HOME/miniforge3/envs/ani_nucdiff/bin/fastANI" \
  --label "fastANI + nucdiff on 5 (reference, hybrid assembly) pairs -> analysis/assembly_comparison/" \
  -- "$HOME/miniforge3/envs/seqqc/bin/python3" tools/compare_assemblies.py
