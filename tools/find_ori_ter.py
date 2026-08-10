#!/usr/bin/env python3
"""Independently locate oriC and terminus on a bacterial chromosome via cumulative
GC-skew (Lobry's method): skew_i = (G-C)/(G+C) per window, cumulative-summed across
the genome. The global minimum of the cumulative curve is oriC, the global maximum
is ter (leading-strand convention). Cross-checked against a user-supplied
approximate oriC region to catch a flipped sign convention.

Usage:
    python find_ori_ter.py <reference.fasta> <output_dir> [--window 1000] [--approx-oric BP]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _genome_utils import read_single_fasta


def compute_cumulative_skew(seq, window):
    n = len(seq)
    positions = []
    cumulative = []
    running = 0.0
    for start in range(0, n, window):
        chunk = seq[start:start + window]
        g = chunk.count("G")
        c = chunk.count("C")
        skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
        running += skew
        positions.append(start + len(chunk) // 2 + 1)  # 1-based midpoint
        cumulative.append(running)
    return positions, cumulative


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("output_dir")
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--approx-oric", type=int, default=None,
                     help="Approximate expected oriC position, used only to check orientation")
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    header, seq = read_single_fasta(args.reference)
    genome_length = len(seq)
    positions, cumulative = compute_cumulative_skew(seq, args.window)

    min_idx = min(range(len(cumulative)), key=lambda i: cumulative[i])
    max_idx = max(range(len(cumulative)), key=lambda i: cumulative[i])
    oric_candidate = positions[min_idx]
    ter_candidate = positions[max_idx]

    # sanity check orientation against the user's approximate oriC, if given
    swapped = False
    if args.approx_oric is not None:
        def circ_dist(a, b):
            d = abs(a - b) % genome_length
            return min(d, genome_length - d)
        if circ_dist(ter_candidate, args.approx_oric) < circ_dist(oric_candidate, args.approx_oric):
            oric_candidate, ter_candidate = ter_candidate, oric_candidate
            swapped = True

    result = {
        "reference": str(args.reference),
        "header": header,
        "genome_length": genome_length,
        "window": args.window,
        "oriC": oric_candidate,
        "ter": ter_candidate,
        "orientation_swapped_vs_approx": swapped,
        "approx_oric_given": args.approx_oric,
    }

    with open(out_dir / "ori_ter.json", "w") as f:
        json.dump(result, f, indent=2)

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(positions, cumulative, color="steelblue", linewidth=1)
    ax.axvline(oric_candidate, color="green", linestyle="--", label=f"oriC ({oric_candidate:,} bp)")
    ax.axvline(ter_candidate, color="red", linestyle="--", label=f"ter ({ter_candidate:,} bp)")
    ax.set_xlabel("Genomic position (bp)")
    ax.set_ylabel("Cumulative GC skew")
    ax.set_title(f"Cumulative GC-skew — {header}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "gc_skew_diagnostic.png", dpi=200)

    print(f"Genome length: {genome_length:,} bp")
    print(f"oriC (cumulative-skew minimum): {oric_candidate:,} bp")
    print(f"ter  (cumulative-skew maximum): {ter_candidate:,} bp")
    if args.approx_oric is not None:
        print(f"Orientation swapped vs. user's approximate oriC ({args.approx_oric:,}): {swapped}")
    print(f"Saved: {out_dir / 'ori_ter.json'}")
    print(f"Saved: {out_dir / 'gc_skew_diagnostic.png'}")


if __name__ == "__main__":
    main()
