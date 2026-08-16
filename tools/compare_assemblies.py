#!/usr/bin/env python3
"""Compare each strain's Unicycler hybrid assembly against its public reference
genome: Average Nucleotide Identity (FastANI) plus SNPs and larger structural/
chromosomal rearrangements (NucDiff).

Usage: python3 compare_assemblies.py

Requires `fastANI` and `nucdiff` on PATH (see tools/compare_assemblies.sh,
which activates the local `ani_nucdiff` conda env before calling this).

Intended to be invoked through scripts/run_logged.py (see
tools/compare_assemblies.sh) so the whole batch is one log entry.
"""
import os
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# fastANI/nucdiff live in the separate `ani_nucdiff` conda env, not the seqqc
# env this script itself runs under (seqqc has python-docx, needed by
# run_logged.py) -- resolve their full binary paths rather than relying on
# PATH/conda activate. Override with ANI_NUCDIFF_BIN if the env lives elsewhere.
ANI_NUCDIFF_BIN = Path(os.environ.get(
    "ANI_NUCDIFF_BIN",
    str(Path.home() / "miniforge3" / "envs" / "ani_nucdiff" / "bin"),
))
FASTANI = str(ANI_NUCDIFF_BIN / "fastANI")
NUCDIFF = str(ANI_NUCDIFF_BIN / "nucdiff")

# nucdiff shells out to nucmer/delta-filter/show-snps/show-coords (from the
# mummer package, installed alongside it in ani_nucdiff) by bare name, so
# that bin/ dir must be on PATH for the subprocess, not just resolvable itself.
NUCDIFF_ENV = dict(os.environ, PATH=f"{ANI_NUCDIFF_BIN}:{os.environ.get('PATH', '')}")

# PairID | reference fasta (relative to PROJECT_ROOT) | hybrid-assembly fasta (relative to PROJECT_ROOT)
# Source filenames are inconsistently prefixed vs. the strain IDs used
# elsewhere in the project (8VU.fna, Hgn_AB042_n.fasta) -- paired here by
# strain identity, not filename.
PAIRS = [
    ("mff", "reference_seq/mff.fna", "hybrid_assembly/Hgn_mff_n.fasta"),
    ("VU", "reference_seq/8VU.fna", "hybrid_assembly/Hgn_VU_n.fasta"),
    ("AB30", "reference_seq/AB30.fna", "hybrid_assembly/Hgn_AB30_n.fasta"),
    ("AB42", "reference_seq/AB42.fna", "hybrid_assembly/Hgn_AB042_n.fasta"),
    ("Lac4", "reference_seq/Lac4.fna", "hybrid_assembly/Hgn_Lac4_n.fasta"),
]

SUMMARY_HEADER = "PairID\tReference\tQuery\tANI\tBidirectionalFragmentMappings\tTotalQueryFragments\tFastANIVersion"


def tool_version(tool, flag="--version"):
    r = subprocess.run([tool, flag], capture_output=True, text=True)
    out = (r.stdout or r.stderr).strip()
    return out or "unknown"


def run_fastani(reference, query, out_tsv, threads):
    r = subprocess.run(
        [FASTANI, "-q", str(query), "-r", str(reference), "-o", str(out_tsv), "--threads", str(threads)],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        return None, f"fastANI failed: {r.stderr[-500:]}"
    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        # genomes too divergent for FastANI to estimate ANI (below ~80% identity)
        return ("NA", "NA", "NA"), None
    line = out_tsv.read_text().splitlines()[0]
    fields = line.split("\t")
    _q, _r, ani, bidir, total = fields[:5]
    return (ani, bidir, total), None


def run_nucdiff(reference, query, out_dir, pair_id):
    out_dir.mkdir(parents=True, exist_ok=True)
    r = subprocess.run(
        [NUCDIFF, str(reference), str(query), str(out_dir), pair_id],
        capture_output=True, text=True, env=NUCDIFF_ENV,
    )
    if r.returncode != 0:
        return False, f"nucdiff failed: {r.stderr[-500:]}"
    return True, None


def process_pair(pair_id, ref_rel, query_rel, out_root, threads):
    """Returns (ani_row_or_None, ani_err_or_None, nucdiff_err_or_None).
    FastANI and NucDiff are independent -- a NucDiff failure must not discard
    an already-successful FastANI result."""
    reference = PROJECT_ROOT / ref_rel
    query = PROJECT_ROOT / query_rel
    missing = [str(p) for p in (reference, query) if not p.is_file()]
    if missing:
        err = f"missing input file(s): {', '.join(missing)}"
        return None, err, err

    out_dir = out_root / pair_id
    out_dir.mkdir(parents=True, exist_ok=True)

    ani_tsv = out_dir / f"{pair_id}_fastani.tsv"
    ani_result, ani_err = run_fastani(reference, query, ani_tsv, threads)
    ani_row = (ref_rel, query_rel) + ani_result if ani_result else None

    nucdiff_dir = out_dir / "nucdiff"
    ok, nucdiff_err = run_nucdiff(reference, query, nucdiff_dir, pair_id)

    return ani_row, ani_err, nucdiff_err


def main():
    threads = os.cpu_count() or 4

    out_root = PROJECT_ROOT / "analysis" / "assembly_comparison"
    out_root.mkdir(parents=True, exist_ok=True)
    summary_tsv = out_root / "summary_ani.tsv"

    fastani_version = tool_version(FASTANI)
    nucdiff_version = tool_version(NUCDIFF)
    print(f"fastANI version in use: {fastani_version}")
    print(f"nucdiff version in use: {nucdiff_version}")
    print(f"Threads: {threads}")
    print()

    rows = [SUMMARY_HEADER]
    ani_failures = []
    nucdiff_failures = []
    for i, (pair_id, ref_rel, query_rel) in enumerate(PAIRS, 1):
        print(f"[{i}/{len(PAIRS)}] {pair_id} ...", flush=True)
        ani_row, ani_err, nucdiff_err = process_pair(pair_id, ref_rel, query_rel, out_root, threads)

        if ani_row:
            _ref, _query, ani, bidir, total = ani_row
            rows.append(f"{pair_id}\t{ref_rel}\t{query_rel}\t{ani}\t{bidir}\t{total}\t{fastani_version}")
            print(f"  ANI={ani}% ({bidir}/{total} fragments mapped)")
        else:
            ani_failures.append((pair_id, ani_err))
            print(f"  FastANI FAILED: {ani_err}")

        if nucdiff_err:
            nucdiff_failures.append((pair_id, nucdiff_err))
            print(f"  NucDiff FAILED: {nucdiff_err}")
        else:
            print(f"  NucDiff OK")

    summary_tsv.write_text("\n".join(rows) + "\n")

    print()
    print(f"=== Assembly comparison summary ===")
    print(f"FastANI: {len(PAIRS) - len(ani_failures)}/{len(PAIRS)} succeeded" + (f", {len(ani_failures)} FAILED" if ani_failures else ""))
    print(f"NucDiff: {len(PAIRS) - len(nucdiff_failures)}/{len(PAIRS)} succeeded" + (f", {len(nucdiff_failures)} FAILED" if nucdiff_failures else ""))
    print(f"ANI summary table: {summary_tsv}")
    print(f"NucDiff outputs in: {out_root}/<PairID>/nucdiff/")
    for pair_id, err in ani_failures:
        print(f"  {pair_id}: FastANI FAILED. {err}")
    for pair_id, err in nucdiff_failures:
        print(f"  {pair_id}: NucDiff FAILED. {err}")

    sys.exit(1 if (ani_failures or nucdiff_failures) else 0)


if __name__ == "__main__":
    main()
