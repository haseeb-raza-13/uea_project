#!/usr/bin/env python3
"""Run a command and automatically append a record of it to docs/QC_Project_Log.docx.

Usage:
    python run_logged.py [--purpose "why this was run"] [--label "text for the log's Command column"] -- <command> [args...]

Example:
    python run_logged.py --purpose "QC of raw short reads" -- \\
        fastqc short_read_data/sample_R1.fastq.gz -o qc_output

--label overrides what's stored in the docx "Command" column (e.g. for batch
commands with dozens of file arguments, where the literal command line would be
unreadable) — the real command still runs exactly as given after '--'.

--tool overrides which binary's version is looked up for the "Tool (version)"
column (e.g. when the executed command is a Python driver script that internally
runs bbduk.sh many times — the column should show BBDuk's version, not Python's).

Note: the command runs directly (no shell), so shell redirection like '>' is
not supported — use each tool's own -o/--output flag instead.
"""
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import _doclog as dl

VERSION_FLAGS = {
    "fastqc": ["--version"],
    "fastp": ["--version"],
    "trimmomatic": ["-version"],
    "cutadapt": ["--version"],
    "bbduk.sh": ["--version"],
    "bbmerge.sh": ["--version"],
    "multiqc": ["--version"],
    "conda": ["--version"],
    "mamba": ["--version"],
}


def get_tool_version(tool):
    name = Path(tool).name
    flags = VERSION_FLAGS.get(name, ["--version"])
    try:
        result = subprocess.run([tool] + flags, capture_output=True, text=True, timeout=30)
        lines = [l for l in (result.stdout + result.stderr).splitlines() if l.strip()]
        if not lines:
            return "(no version output)"
        # skip launcher/echo lines (e.g. bbduk.sh prints its underlying `java -cp ... --version` invocation)
        candidates = [l for l in lines if not l.lstrip().startswith("java ") and "-cp " not in l]
        pool = candidates or lines
        version_line = next((l for l in pool if re.search(r"\bversion\b", l, re.IGNORECASE)), pool[0])
        return version_line.strip()[:150]
    except Exception as exc:
        return f"(version check failed: {exc})"


def parse_args(argv):
    if "--" not in argv:
        print("Usage: run_logged.py [--purpose TEXT] -- <command> [args...]", file=sys.stderr)
        sys.exit(2)
    sep = argv.index("--")
    pre, cmd = argv[:sep], argv[sep + 1:]
    if not cmd:
        print("No command given after --", file=sys.stderr)
        sys.exit(2)
    purpose = ""
    label = ""
    tool = ""
    i = 0
    while i < len(pre):
        if pre[i] == "--purpose" and i + 1 < len(pre):
            purpose = pre[i + 1]
            i += 2
        elif pre[i] == "--label" and i + 1 < len(pre):
            label = pre[i + 1]
            i += 2
        elif pre[i] == "--tool" and i + 1 < len(pre):
            tool = pre[i + 1]
            i += 2
        else:
            i += 1
    return purpose, label, tool, cmd


def main():
    purpose, label, tool, cmd = parse_args(sys.argv[1:])
    tool_version = get_tool_version(tool or cmd[0])

    start = datetime.now()
    t0 = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        exit_code = result.returncode
        output = (result.stdout or "") + (result.stderr or "")
    except FileNotFoundError as exc:
        exit_code = 127
        output = str(exc)
    duration = time.time() - t0

    status = "Success" if exit_code == 0 else f"Failed (exit {exit_code})"
    raw_lines = output.strip().splitlines()
    # some tools (e.g. FastQC in multi-threaded mode) spam stray "null" lines to stderr
    lines = [l for l in raw_lines if l.strip() and l.strip().lower() != "null"]
    omitted = len(raw_lines) - len(lines)
    if len(lines) > 25:
        header = f"... [showing last 25 of {len(lines)} output line(s)"
        header += f", {omitted} filler line(s) omitted" if omitted else ""
        header += "] ..."
        notes = header + "\n" + "\n".join(lines[-25:])
    else:
        notes = "\n".join(lines)
    notes = (notes + f"\n(duration: {duration:.1f}s)").strip()

    doc = dl.load_doc()
    step = dl.append_row(
        doc,
        timestamp=start,
        purpose=purpose or "(not specified)",
        tool_version=tool_version,
        command=label or " ".join(shlex.quote(c) for c in cmd),
        status=status,
        notes=notes,
    )

    print(output)
    print(f"[logged as step {step} in {dl.DOC_PATH.name}]", file=sys.stderr)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
