#!/usr/bin/env python3
"""Parse Sample Name illumina.xlsx into a tidy per-sample metadata CSV
(strain, replicate, treatment, growth phase) keyed by the PID-2861-<N> sample id
used throughout the pipeline's file naming. Project-wide (not mff-specific) so it
can be reused for the other 4 strains later.
"""
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
XLSX_PATH = PROJECT_ROOT / "short_read_data" / "short_read_seq_1" / "Sample Name illumina.xlsx"
OUT_PATH = PROJECT_ROOT / "mapping" / "coverage" / "sample_metadata.csv"


def main():
    df = pd.read_excel(XLSX_PATH)
    rows = []
    for _, row in df.iterrows():
        tokens = str(row["Sample ID"]).split()
        if len(tokens) != 4:
            print(f"WARNING: unexpected Sample ID format, skipping: {row['Sample ID']!r}", file=sys.stderr)
            continue
        strain, replicate, treatment, phase = tokens
        sample_number = int(row["Sample number"])
        rows.append({
            "sample_number": sample_number,
            "pid_sample": f"PID-2861-{sample_number}",
            "sample_id_raw": row["Sample ID"],
            "strain": strain,
            "replicate": replicate,
            "treatment": treatment,
            "phase": phase,
        })

    out_df = pd.DataFrame(rows).sort_values("sample_number")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(out_df)} sample records to {OUT_PATH}")


if __name__ == "__main__":
    main()
