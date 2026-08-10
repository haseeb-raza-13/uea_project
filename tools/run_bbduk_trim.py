#!/usr/bin/env python3
"""Remove Illumina/Nextera adapter contamination from all short-read samples using BBDuk.

Pairs R1/R2 FASTQ files by sample name across short_read_data/short_read_seq_1/2/3
(mates are scattered arbitrarily across those three folders), then runs BBDuk's
standard adapter-trimming preset on each pair. Intended to be invoked through
scripts/run_logged.py (see tools/run_bbduk_trim.sh) so the whole batch is recorded
as one log entry, not one per sample.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SHORT_READ_DIRS = [
    PROJECT_ROOT / "short_read_data" / "short_read_seq_1",
    PROJECT_ROOT / "short_read_data" / "short_read_seq_2",
    PROJECT_ROOT / "short_read_data" / "short_read_seq_3",
]
OUT_DIR = PROJECT_ROOT / "trimmed_data" / "short_read"
STATS_DIR = OUT_DIR / "bbduk_stats"

FASTQ_RE = re.compile(r"^(?P<sample>.+)_S\d+_R(?P<read>[12])_001\.fastq(\.gz)?$")

SUMMARY_LINE_RE = re.compile(
    r"^(?P<label>Input|KTrimmed|Trimmed by overlap|Total Removed|Result):\s+"
    r"(?P<reads>[\d,]+)\s+reads\s*(?:\([\d.]+%\))?\s+"
    r"(?P<bases>[\d,]+)\s+bases\s*(?:\([\d.]+%\))?"
)


def find_pairs():
    files_by_read = {"1": {}, "2": {}}
    for d in SHORT_READ_DIRS:
        for f in d.glob("*.fastq*"):
            m = FASTQ_RE.match(f.name)
            if not m:
                continue
            sample, read = m.group("sample"), m.group("read")
            files_by_read[read][sample] = f

    samples = sorted(set(files_by_read["1"]) | set(files_by_read["2"]))
    pairs = {}
    incomplete = []
    for sample in samples:
        r1 = files_by_read["1"].get(sample)
        r2 = files_by_read["2"].get(sample)
        if r1 and r2:
            pairs[sample] = (r1, r2)
        else:
            incomplete.append((sample, bool(r1), bool(r2)))

    if incomplete:
        print("ERROR: incomplete R1/R2 pairs found:", file=sys.stderr)
        for sample, has_r1, has_r2 in incomplete:
            print(f"  {sample}: R1={'yes' if has_r1 else 'MISSING'} R2={'yes' if has_r2 else 'MISSING'}", file=sys.stderr)
        sys.exit(1)

    return pairs


def parse_bbduk_summary(text):
    totals = {}
    for line in text.splitlines():
        m = SUMMARY_LINE_RE.match(line.strip())
        if m:
            totals[m.group("label")] = (
                int(m.group("reads").replace(",", "")),
                int(m.group("bases").replace(",", "")),
            )
    return totals


def main():
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        print("ERROR: CONDA_PREFIX not set — run this inside the 'seqqc' conda environment.", file=sys.stderr)
        sys.exit(1)
    adapters_ref = Path(conda_prefix) / "share" / "bbmap" / "resources" / "adapters.fa"
    if not adapters_ref.exists():
        print(f"ERROR: adapter reference not found at {adapters_ref}", file=sys.stderr)
        sys.exit(1)

    threads = max(os.cpu_count() - 2, 1) if os.cpu_count() else 4

    pairs = find_pairs()
    print(f"Found {len(pairs)} complete R1/R2 sample pairs.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    STATS_DIR.mkdir(parents=True, exist_ok=True)

    failures = []
    total_input_bases = 0
    total_removed_bases = 0

    for i, (sample, (r1, r2)) in enumerate(sorted(pairs.items()), 1):
        out1 = OUT_DIR / f"{sample}_R1.trimmed.fastq.gz"
        out2 = OUT_DIR / f"{sample}_R2.trimmed.fastq.gz"
        stats_file = STATS_DIR / f"{sample}_stats.txt"

        cmd = [
            "bbduk.sh",
            f"in={r1}",
            f"in2={r2}",
            f"out={out1}",
            f"out2={out2}",
            f"ref={adapters_ref}",
            "ktrim=r", "k=23", "mink=11", "hdist=1", "tpe", "tbo",
            f"threads={threads}",
            f"stats={stats_file}",
        ]
        print(f"[{i}/{len(pairs)}] {sample} ...", flush=True)
        result = subprocess.run(cmd, capture_output=True, text=True)
        combined = result.stdout + result.stderr

        if result.returncode != 0:
            failures.append(sample)
            print(f"  FAILED (exit {result.returncode})")
            continue

        summary = parse_bbduk_summary(combined)
        if "Input" in summary:
            total_input_bases += summary["Input"][1]
        removed = summary.get("Total Removed")
        if removed is None:
            kt = summary.get("KTrimmed", (0, 0))
            ov = summary.get("Trimmed by overlap", (0, 0))
            removed = (kt[0] + ov[0], kt[1] + ov[1])
        total_removed_bases += removed[1]

    ok = len(pairs) - len(failures)
    pct = (total_removed_bases / total_input_bases * 100) if total_input_bases else 0.0
    print()
    print("=== BBDuk adapter-trimming summary ===")
    print(f"Samples processed: {ok}/{len(pairs)} succeeded" + (f", {len(failures)} FAILED: {failures}" if failures else ""))
    print(f"Total input bases: {total_input_bases:,}")
    print(f"Total bases removed (adapter/overlap trimming): {total_removed_bases:,} ({pct:.2f}%)")
    print(f"Trimmed reads written to: {OUT_DIR}")
    print(f"Per-sample BBDuk stats written to: {STATS_DIR}")

    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
