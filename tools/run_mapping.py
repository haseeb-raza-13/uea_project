#!/usr/bin/env python3
"""Align a strain's trimmed sample pairs to its reference (Bowtie2), then
deduplicate (samtools sort/fixmate/markdup). PCR duplicates must be removed before
marker-frequency analysis, since they create false coverage spikes unrelated to
real DNA replication.

Usage: python3 run_mapping.py <strain>

The sample list for the strain is read from mapping/coverage/sample_metadata.csv
(built by tools/build_sample_metadata.py).

Intended to be invoked through scripts/run_logged.py (see tools/run_mapping.sh) so
the whole per-strain batch is recorded as one log entry.
"""
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TRIMMED_DIR = PROJECT_ROOT / "trimmed_data" / "short_read"
METADATA_CSV = PROJECT_ROOT / "mapping" / "coverage" / "sample_metadata.csv"

ALIGN_RATE_RE = re.compile(r"([\d.]+)% overall alignment rate")


def samples_for_strain(strain):
    with open(METADATA_CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    return sorted(
        (r["pid_sample"] for r in rows if r["strain"] == strain),
        key=lambda s: int(s.rsplit("-", 1)[1]),
    )


def run(cmd, **kwargs):
    return subprocess.run(cmd, capture_output=True, text=True, **kwargs)


def parse_markdup_stats(text):
    stats = {}
    for line in text.splitlines():
        m = re.match(r"^([A-Z][A-Z_ ]*[A-Z]):\s*(\d+)", line.strip())
        if m:
            stats[m.group(1)] = int(m.group(2))
    return stats


def process_sample(sample, index_prefix, bam_dir, threads):
    r1 = TRIMMED_DIR / f"{sample}_R1.trimmed.fastq.gz"
    r2 = TRIMMED_DIR / f"{sample}_R2.trimmed.fastq.gz"
    if not r1.exists() or not r2.exists():
        return {"sample": sample, "ok": False, "error": f"missing trimmed fastq(s) for {sample}"}

    sam_path = bam_dir / f"{sample}.sam"
    sorted_bam = bam_dir / f"{sample}.sorted.bam"
    namesorted_bam = bam_dir / f"{sample}.namesorted.bam"
    fixmate_bam = bam_dir / f"{sample}.fixmate.bam"
    coordsorted_bam = bam_dir / f"{sample}.fixmate.sorted.bam"
    dedup_bam = bam_dir / f"{sample}.dedup.bam"

    bt2 = run([
        "bowtie2", "-x", str(index_prefix), "-1", str(r1), "-2", str(r2),
        "-p", str(threads), "--rg-id", sample, "--rg", f"SM:{sample}",
        "-S", str(sam_path),
    ])
    if bt2.returncode != 0:
        return {"sample": sample, "ok": False, "error": f"bowtie2 failed: {bt2.stderr[-500:]}"}
    m = ALIGN_RATE_RE.search(bt2.stderr)
    align_rate = float(m.group(1)) if m else None

    r = run(["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), str(sam_path)])
    sam_path.unlink(missing_ok=True)
    if r.returncode != 0:
        return {"sample": sample, "ok": False, "error": f"samtools sort failed: {r.stderr[-500:]}"}

    r = run(["samtools", "sort", "-n", "-@", str(threads), "-o", str(namesorted_bam), str(sorted_bam)])
    if r.returncode != 0:
        return {"sample": sample, "ok": False, "error": f"samtools sort -n failed: {r.stderr[-500:]}"}

    r = run(["samtools", "fixmate", "-m", str(namesorted_bam), str(fixmate_bam)])
    namesorted_bam.unlink(missing_ok=True)
    if r.returncode != 0:
        return {"sample": sample, "ok": False, "error": f"samtools fixmate failed: {r.stderr[-500:]}"}

    r = run(["samtools", "sort", "-@", str(threads), "-o", str(coordsorted_bam), str(fixmate_bam)])
    fixmate_bam.unlink(missing_ok=True)
    if r.returncode != 0:
        return {"sample": sample, "ok": False, "error": f"samtools sort (post-fixmate) failed: {r.stderr[-500:]}"}

    r = run(["samtools", "markdup", "-r", "-s", str(coordsorted_bam), str(dedup_bam)])
    coordsorted_bam.unlink(missing_ok=True)
    if r.returncode != 0:
        return {"sample": sample, "ok": False, "error": f"samtools markdup failed: {r.stderr[-500:]}"}
    dup_stats = parse_markdup_stats(r.stderr)

    r = run(["samtools", "index", str(dedup_bam)])
    if r.returncode != 0:
        return {"sample": sample, "ok": False, "error": f"samtools index failed: {r.stderr[-500:]}"}

    sorted_bam.unlink(missing_ok=True)

    examined = dup_stats.get("EXAMINED", 0)
    dup_total = dup_stats.get("DUPLICATE TOTAL", 0)
    dup_pct = (dup_total / examined * 100) if examined else 0.0

    return {
        "sample": sample, "ok": True,
        "align_rate": align_rate,
        "examined": examined,
        "dup_total": dup_total,
        "dup_pct": dup_pct,
        "bam": str(dedup_bam),
    }


def main():
    if len(sys.argv) != 2:
        print("Usage: run_mapping.py <strain>", file=sys.stderr)
        sys.exit(2)
    strain = sys.argv[1]

    index_prefix = PROJECT_ROOT / "mapping" / "reference" / strain / f"{strain}_index"
    bam_dir = PROJECT_ROOT / "mapping" / "bam" / strain

    if not Path(str(index_prefix) + ".1.bt2").exists():
        print(f"ERROR: Bowtie2 index not found at {index_prefix}. Run tools/build_reference_index.sh first.", file=sys.stderr)
        sys.exit(1)

    samples = samples_for_strain(strain)
    if not samples:
        print(f"ERROR: no samples found for strain '{strain}' in {METADATA_CSV}", file=sys.stderr)
        sys.exit(1)

    threads = max(os.cpu_count() - 2, 1) if os.cpu_count() else 4
    bam_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for i, sample in enumerate(samples, 1):
        print(f"[{i}/{len(samples)}] {sample} ...", flush=True)
        res = process_sample(sample, index_prefix, bam_dir, threads)
        results.append(res)
        if not res["ok"]:
            print(f"  FAILED: {res['error']}")
        else:
            print(f"  align rate={res['align_rate']}%  duplicates={res['dup_pct']:.2f}% ({res['dup_total']}/{res['examined']})")

    failures = [r for r in results if not r["ok"]]
    ok = [r for r in results if r["ok"]]

    print()
    print(f"=== Mapping + dedup summary ({strain}, {len(samples)} samples) ===")
    fail_suffix = f", {len(failures)} FAILED" if failures else ""
    print(f"Samples processed: {len(ok)}/{len(samples)} succeeded{fail_suffix}")
    for r in ok:
        print(f"  {r['sample']}: alignment {r['align_rate']}%, duplicates removed {r['dup_pct']:.2f}%")
    for r in failures:
        print(f"  {r['sample']}: FAILED. {r['error']}")
    if ok:
        avg_align = sum(r["align_rate"] for r in ok if r["align_rate"] is not None) / len(ok)
        avg_dup = sum(r["dup_pct"] for r in ok) / len(ok)
        print(f"Average alignment rate: {avg_align:.2f}%")
        print(f"Average duplicate rate removed: {avg_dup:.2f}%")
    print(f"Deduplicated BAMs written to: {bam_dir}")

    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
