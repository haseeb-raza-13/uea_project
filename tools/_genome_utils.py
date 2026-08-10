"""Shared circular-genome coordinate helpers for the mapping/replication-profile pipeline.

Not meant to be run directly.
"""
from pathlib import Path


def read_single_fasta(path):
    """Read a single-record FASTA, return (header, sequence) with sequence upper-cased."""
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"{path} contains more than one FASTA record; expected exactly one chromosome.")
                header = line[1:].strip()
            else:
                seq_parts.append(line.strip())
    if header is None:
        raise ValueError(f"No FASTA records found in {path}")
    return header, "".join(seq_parts).upper()


def replichore_normalize(pos, oric, ter, genome_length):
    """Map a 1-based genomic position onto [-1, +1], oriC=0, ter=+1 on the forward
    replichore and -1 on the reverse replichore, following Skovgaard et al. 2011.
    """
    fwd_total = (ter - oric) % genome_length
    bwd_total = genome_length - fwd_total
    fwd = (pos - oric) % genome_length
    if fwd_total == 0 or bwd_total == 0:
        return 0.0
    if fwd <= fwd_total:
        return fwd / fwd_total
    else:
        bwd = genome_length - fwd
        return -bwd / bwd_total
