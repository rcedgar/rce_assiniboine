#!/usr/bin/env python3
"""Write per-sequence ACGT frequencies for a DNA FASTA file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import fasta

_BASES = ("A", "C", "G", "T")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Write DNA base frequencies for each FASTA record."
    )
    parser.add_argument("fasta", help="Input FASTA file")
    args = parser.parse_args(argv)

    def onseq(label: str, seq: str) -> None:
        seq_u = seq.upper()
        length = len(seq_u)
        counts = [seq_u.count(base) for base in _BASES]
        fields = [label, str(length)]
        fields.extend("%.3f" % (count / length) for count in counts)
        sys.stdout.write("\t".join(fields) + "\n")

    fasta.ReadSeqsOnSeq(args.fasta, onseq)


if __name__ == "__main__":
    main()
