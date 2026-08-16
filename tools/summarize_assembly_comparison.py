#!/usr/bin/env python3
"""Consolidate per-strain FastANI + NucDiff assembly-comparison output into one
cross-strain summary table.

Usage: python3 summarize_assembly_comparison.py

Reads analysis/assembly_comparison/summary_ani.tsv (written by
tools/compare_assemblies.sh) for ANI/fragment-mapping stats, and each
strain's analysis/assembly_comparison/<strain>/nucdiff/results/<strain>_stat.out
for NucDiff's SNP and structural-rearrangement counts. The stat.out format
(plain "Key<TAB>Value" lines, blank-line-separated sections, verified against
real output for all 5 strains) is:

    Total number	58
    Insertions	9
    Deletions	1
    Substitutions	43
    Translocations	0
    Relocations	2
    Reshufflings	1
    Reshuffled blocks	2
    Inversions	0
    Unaligned sequences	2

    Uncovered ref regions num	0
    ...
    DETAILED INFORMATION:
    ...

Only the top summary block (before the blank line) is parsed here.

Intended to be invoked through scripts/run_logged.py (see
tools/summarize_assembly_comparison.sh).
"""
import csv
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STRAINS = ["mff", "VU", "AB30", "AB42", "Lac4"]

STAT_FIELDS = [
    ("Total number", "total_differences"),
    ("Insertions", "insertions"),
    ("Deletions", "deletions"),
    ("Substitutions", "substitutions"),
    ("Translocations", "translocations"),
    ("Relocations", "relocations"),
    ("Reshufflings", "reshufflings"),
    ("Reshuffled blocks", "reshuffled_blocks"),
    ("Inversions", "inversions"),
    ("Unaligned sequences", "unaligned_sequences"),
]


def load_ani_summary(ani_tsv):
    rows = {}
    with open(ani_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows[row["PairID"]] = row
    return rows


def contig_lengths(fasta_path):
    lengths = []
    cur = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lengths.append(cur)
    return lengths


def assembly_completeness(query_rel_path):
    lengths = contig_lengths(PROJECT_ROOT / query_rel_path)
    total = sum(lengths)
    n_contigs = len(lengths)
    largest = max(lengths) if lengths else 0
    frac = largest / total if total else float("nan")
    return {
        "n_contigs": n_contigs,
        "total_assembly_len": total,
        "largest_contig_len": largest,
        "largest_contig_frac": f"{frac:.4f}" if total else "NA",
    }


def load_nucdiff_stats(strain, out_root):
    stat_path = out_root / strain / "nucdiff" / "results" / f"{strain}_stat.out"
    if not stat_path.exists():
        return None, f"{stat_path} not found"

    parsed = {}
    for line in stat_path.read_text().splitlines():
        if not line.strip():
            break  # end of the top summary block
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        parsed[key.strip()] = value.strip()

    stats = {}
    missing_keys = []
    for label, col in STAT_FIELDS:
        if label not in parsed:
            missing_keys.append(label)
            stats[col] = "NA"
        else:
            stats[col] = parsed[label]

    err = f"missing expected key(s) in {stat_path.name}: {', '.join(missing_keys)}" if missing_keys else None
    return stats, err


def main():
    out_root = PROJECT_ROOT / "analysis" / "assembly_comparison"
    ani_tsv = out_root / "summary_ani.tsv"

    if not ani_tsv.exists():
        print(f"ERROR: {ani_tsv} not found. Run tools/compare_assemblies.sh first.",
              file=sys.stderr)
        sys.exit(1)

    ani_rows = load_ani_summary(ani_tsv)

    out_csv = out_root / "summary_table.csv"
    completeness_fields = ["n_contigs", "total_assembly_len", "largest_contig_len", "largest_contig_frac"]
    fieldnames = [
        "strain", "reference", "query", "ani_percent",
        "bidirectional_fragments", "total_query_fragments", "fraction_mapped",
    ] + [col for _, col in STAT_FIELDS] + completeness_fields

    missing_strains = []
    warnings = []
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for strain in STRAINS:
            if strain not in ani_rows:
                missing_strains.append(strain)
                continue
            r = ani_rows[strain]
            ani = r["ANI"]
            bidir = r["BidirectionalFragmentMappings"]
            total = r["TotalQueryFragments"]
            try:
                fraction_mapped = f"{int(bidir) / int(total):.4f}"
            except (ValueError, ZeroDivisionError):
                fraction_mapped = "NA"

            stats, err = load_nucdiff_stats(strain, out_root)
            if err:
                warnings.append(f"{strain}: {err}")
            if stats is None:
                stats = {col: "NA" for _, col in STAT_FIELDS}

            row = {
                "strain": strain,
                "reference": r["Reference"],
                "query": r["Query"],
                "ani_percent": ani,
                "bidirectional_fragments": bidir,
                "total_query_fragments": total,
                "fraction_mapped": fraction_mapped,
            }
            row.update(stats)
            row.update(assembly_completeness(r["Query"]))
            writer.writerow(row)

    print(f"=== Assembly comparison summary ===")
    print(f"Strains summarized: {len(STRAINS) - len(missing_strains)}/{len(STRAINS)}")
    if missing_strains:
        print(f"  Missing from {ani_tsv.name}: {', '.join(missing_strains)}")
    if warnings:
        print(f"NucDiff parsing warnings:")
        for w in warnings:
            print(f"  {w}")
    print(f"Table written to: {out_csv}")

    sys.exit(1 if missing_strains else 0)


if __name__ == "__main__":
    main()
