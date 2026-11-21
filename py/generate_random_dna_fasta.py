#!/usr/bin/env python3
"""Write one random DNA FASTA record with prescribed base frequencies."""

from __future__ import annotations

import argparse
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import fasta

_BASES = ("A", "C", "G", "T")


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def base_counts(length: int, freqs: list[float]) -> list[int]:
    """Integer ACGT counts summing to length, each within < 1 of its quota.

    Frequencies are normalized to sum to 1, then split by the largest-remainder
    method: floor(p * L), with leftover bases given to the largest fractional parts.
    """
    total = sum(freqs)
    scaled = [freq / total * length for freq in freqs]
    counts = [int(math.floor(value)) for value in scaled]
    leftover = length - sum(counts)
    if leftover >= 0:
        order = sorted(
            range(len(counts)),
            key=lambda i: (scaled[i] - counts[i], i),
            reverse=True,
        )
        if leftover > len(order):
            _die("internal error: cannot place leftover bases within rounding error")
        for index in range(leftover):
            counts[order[index]] += 1
    else:
        order = sorted(
            range(len(counts)),
            key=lambda i: (scaled[i] - counts[i], i),
        )
        for index in order:
            if leftover == 0:
                break
            if counts[index] == 0:
                continue
            counts[index] -= 1
            leftover += 1
        if leftover != 0:
            _die("internal error: cannot remove extra bases within rounding error")
    if sum(counts) != length:
        _die("internal error: base counts do not sum to L")
    for count, value in zip(counts, scaled):
        if count < 0 or abs(count - value) >= 1.0:
            _die("internal error: base count is outside rounding error")
    return counts


def random_sequence(counts: list[int]) -> str:
    bases: list[str] = []
    for base, count in zip(_BASES, counts):
        bases.extend([base] * count)
    random.shuffle(bases)
    return "".join(bases)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Write N random DNA FASTA records with prescribed base frequencies."
    )
    parser.add_argument("label", help="FASTA defline label")
    parser.add_argument("L", type=int, help="Sequence length")
    parser.add_argument("fA", type=float, help="Frequency of A")
    parser.add_argument("fC", type=float, help="Frequency of C")
    parser.add_argument("fG", type=float, help="Frequency of G")
    parser.add_argument("fT", type=float, help="Frequency of T")
    parser.add_argument("seed", type=float, help="Random number seed (integer)")
    args = parser.parse_args(argv)

    if args.L < 1:
        _die("L must be >= 1")

    freqs = [args.fA, args.fC, args.fG, args.fT]
    for base, freq in zip(_BASES, freqs):
        if not math.isfinite(freq) or freq < 0.0:
            _die(f"f{base} must be finite and >= 0")

    total = sum(freqs)
    if total < 0.99 or total > 1.01:
        _die(f"frequencies sum to {total}, expected 0.99..1.01")

    random.seed(int(args.seed))
    seq = random_sequence(base_counts(args.L, freqs))
    fasta.WriteSeq(sys.stdout, seq, args.label)


if __name__ == "__main__":
    main()
