#!/usr/bin/env python3
"""Score a multi-FASTA with circular mononucleotide-shuffle MFE Z-scores.

Each record is folded with RNAfold --circ, together with mononucleotide
shuffles of the same sequence. Z is

    (original MFE - mean) / sample standard deviation

of the shuffled MFEs. Output is label and Z.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import statistics
import sys
from pathlib import Path

import permute_fold


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return (id, sequence) for every FASTA record. Id is the first header token."""
    records: list[tuple[str, str]] = []
    seq_id = ""
    parts: list[str] = []

    def finish() -> None:
        if seq_id or parts:
            records.append((seq_id, "".join(parts)))

    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if seq_id or parts:
                        finish()
                    header = line[1:].strip()
                    seq_id = header.split()[0] if header else ""
                    parts = []
                else:
                    parts.append(line.replace(" ", "").replace("\t", ""))
        if seq_id or parts:
            finish()
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    if not records:
        _die(f"ERROR: no FASTA records in {path}")
    return records


def _prepare(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _label(seq_id: str, index: int) -> str:
    label = seq_id or f"record_{index}"
    return label.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _is_acgt(seq: str) -> bool:
    return bool(seq) and all(base in "ACGT" for base in seq)


def _fmt(value: float | None) -> str:
    if not isinstance(value, float) or not math.isfinite(value):
        return ""
    return f"{value:.6g}"


def _measured_z(seq_id: str, seq: str, n_permutations: int, seed: int, threads: int) -> float | None:
    permuted = permute_fold.permute_sequences(seq, seq_id, n_permutations, seed)
    fasta_text = permute_fold.format_fasta([(seq_id, seq)])
    fasta_text += permute_fold.format_fasta(permuted)
    try:
        fold_text = permute_fold.run_rnafold(fasta_text, threads)
    except SystemExit:
        return None
    _structure, mfe_rows = permute_fold.parse_fold(fold_text)
    if len(mfe_rows) != n_permutations + 1:
        print(
            f"ERROR: expected {n_permutations + 1} MFE values, parsed {len(mfe_rows)}",
            file=sys.stderr,
        )
        return None
    try:
        original_mfe = float(mfe_rows[0][1])
        shuffled = [float(energy) for _perm_label, energy in mfe_rows[1:]]
    except ValueError:
        return None
    if len(shuffled) < 2:
        return None
    return permute_fold.z_score(original_mfe, statistics.mean(shuffled), statistics.stdev(shuffled))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Write RNAfold shuffle Z-scores for a multi-FASTA.")
    parser.add_argument("--input", required=True, type=Path, help="Multi-FASTA")
    parser.add_argument("--output", required=True, type=Path, help="TSV: label Z")
    parser.add_argument("--permutations", type=int, default=100, help="Mononucleotide shuffles per sequence")
    parser.add_argument("--seed", type=int, default=1, help="Shuffle seed")
    parser.add_argument("--threads", type=int, default=64, help="Maximum RNAfold --jobs")
    args = parser.parse_args(argv)
    if args.permutations < 1:
        _die("ERROR: permutations must be a positive integer")
    if args.seed < 0:
        _die("ERROR: seed must be a non-negative integer")
    if args.threads < 1:
        _die("ERROR: threads must be a positive integer")
    if shutil.which("RNAfold") is None:
        _die("ERROR: RNAfold not found in PATH")

    records = _read_fasta(args.input)
    cpu = os.cpu_count() or 1
    threads = max(1, min(args.threads, cpu, args.permutations + 1))
    total = len(records)
    print(
        f"sequences {total} permutations {args.permutations} RNAfold_jobs {threads}",
        file=sys.stderr,
        flush=True,
    )
    scored = 0
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("label\tZ\n")
        handle.flush()
        for index, (seq_id, raw_seq) in enumerate(records, start=1):
            label = _label(seq_id, index)
            seq = _prepare(raw_seq)
            z_value = None
            if _is_acgt(seq):
                z_value = _measured_z(seq_id, seq, args.permutations, args.seed, threads)
            z_text = _fmt(z_value)
            if z_text:
                scored += 1
            handle.write(f"{label}\t{z_text}\n")
            handle.flush()
            print(f"{index}/{total}\t{label}\t{z_text}", file=sys.stderr, flush=True)
    print(f"scored {scored}", file=sys.stderr, flush=True)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
