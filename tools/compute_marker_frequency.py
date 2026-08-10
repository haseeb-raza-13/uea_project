#!/usr/bin/env python3
"""Compute windowed marker-frequency (read-start count) profiles from deduplicated
BAMs, normalize to the genome-wide median, and fold genomic position onto the
oriC=0 / ter=+-1 replichore-normalized coordinate (Skovgaard et al. 2011 methodology).

Uses `bedtools genomecov -5 -bg` to get a read-start depth track (each read
contributes exactly 1 unit of depth at its 5' position), then bins that track into
fixed-size windows — summing depth*overlap_length over a window equals the total
number of read starts in that window.

Intended to be invoked through scripts/run_logged.py (see
tools/compute_marker_frequency.sh) so the whole batch is one log entry.
"""
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _genome_utils import replichore_normalize

BAM_DIR = PROJECT_ROOT / "mapping" / "bam" / "mff"
OUT_DIR = PROJECT_ROOT / "mapping" / "coverage" / "mff"
ORI_TER_JSON = PROJECT_ROOT / "mapping" / "reference" / "mff" / "ori_ter.json"
WINDOW_SIZES = [1000, 10000]


def bin_bedgraph(bedgraph_text, window_size, genome_length):
    num_windows = genome_length // window_size
    counts = [0.0] * num_windows
    for line in bedgraph_text.splitlines():
        if not line.strip():
            continue
        _, start_s, end_s, depth_s = line.split()
        start, end, depth = int(start_s), int(end_s), float(depth_s)
        if depth == 0:
            continue
        w_start = start // window_size
        w_end = (end - 1) // window_size
        pos = start
        for w in range(w_start, w_end + 1):
            boundary = min((w + 1) * window_size, end)
            overlap = boundary - pos
            if w < num_windows:
                counts[w] += depth * overlap
            pos = boundary
    return counts


def process_sample(sample, bam_path, oric, ter, genome_length):
    r = subprocess.run(["bedtools", "genomecov", "-5", "-bg", "-ibam", str(bam_path)],
                        capture_output=True, text=True)
    if r.returncode != 0:
        return None, f"bedtools genomecov failed: {r.stderr[-500:]}"

    for window_size in WINDOW_SIZES:
        counts = bin_bedgraph(r.stdout, window_size, genome_length)
        nonzero = [c for c in counts if c > 0]
        if not nonzero:
            return None, f"no nonzero windows at {window_size}bp for {sample}"
        baseline = statistics.median(nonzero)  # robust to repeat-region spikes (paper excludes repeat windows explicitly; median substitutes for that here)

        rows = []
        for i, count in enumerate(counts):
            if count <= 0:
                continue
            w_start = i * window_size
            w_end = w_start + window_size
            mid = w_start + window_size // 2 + 1  # 1-based midpoint
            norm_pos = replichore_normalize(mid, oric, ter, genome_length)
            log2_ratio = math.log2(count / baseline)
            rows.append((w_start + 1, w_end, mid, count, norm_pos, log2_ratio))

        out_path = OUT_DIR / f"{sample}_{window_size // 1000}kb.tsv"
        with open(out_path, "w") as f:
            f.write("window_start\twindow_end\twindow_mid\traw_count\tnorm_position\tlog2_ratio\n")
            for row in rows:
                f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]:.1f}\t{row[4]:.6f}\t{row[5]:.6f}\n")

    return True, None


def main():
    if not ORI_TER_JSON.exists():
        print(f"ERROR: {ORI_TER_JSON} not found — run tools/find_ori_ter.sh first.", file=sys.stderr)
        sys.exit(1)
    with open(ORI_TER_JSON) as f:
        ori_ter = json.load(f)
    oric, ter, genome_length = ori_ter["oriC"], ori_ter["ter"], ori_ter["genome_length"]

    bams = sorted(BAM_DIR.glob("*.dedup.bam"))
    if not bams:
        print(f"ERROR: no *.dedup.bam files found in {BAM_DIR} — run tools/run_mapping.sh first.", file=sys.stderr)
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"oriC={oric:,}  ter={ter:,}  genome_length={genome_length:,}")
    print(f"Windows: {WINDOW_SIZES}")

    failures = []
    ok_samples = []
    for i, bam in enumerate(bams, 1):
        sample = bam.stem.replace(".dedup", "")
        print(f"[{i}/{len(bams)}] {sample} ...", flush=True)
        ok, err = process_sample(sample, bam, oric, ter, genome_length)
        if ok:
            ok_samples.append(sample)
        else:
            failures.append((sample, err))
            print(f"  FAILED: {err}")

    print()
    print("=== Marker-frequency computation summary (mff) ===")
    print(f"Samples processed: {len(ok_samples)}/{len(bams)} succeeded" + (f", {len(failures)} FAILED" if failures else ""))
    print(f"Window sizes: {WINDOW_SIZES} bp")
    print(f"oriC={oric:,} bp, ter={ter:,} bp (from {ORI_TER_JSON.name})")
    print(f"TSVs written to: {OUT_DIR}")
    for sample, err in failures:
        print(f"  {sample}: FAILED — {err}")

    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
