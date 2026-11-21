#!/usr/bin/env python3
"""Generate 100 random DNA sequences and check reported base frequencies."""

from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import fasta

_BASES = ("A", "C", "G", "T")
_N = 100
_SEED = 1


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )


def _check_fatal(generate_py: Path) -> None:
    proc = _run(
        [
            sys.executable,
            str(generate_py),
            "bad",
            "100",
            "0.5",
            "0.5",
            "0.5",
            "0.0",
            "1",
        ]
    )
    if proc.returncode == 0:
        _die("frequencies summing to 1.5 were accepted")
    if proc.stdout != "":
        _die("fatal frequency error wrote to stdout")


def _draw_specs(rng: random.Random) -> list[tuple[str, int, list[float]]]:
    specs: list[tuple[str, int, list[float]]] = []
    for index in range(_N):
        length = rng.randint(30, 4000)
        weights = [rng.randint(1, 100) for _ in _BASES]
        total = sum(weights)
        freqs = [weight / total for weight in weights]
        specs.append((f"s{index:03d}", length, freqs))
    return specs


def _generate(generate_py: Path, specs: list[tuple[str, int, list[float]]], fasta_path: Path) -> None:
    with fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
        for label, length, freqs in specs:
            cmd = [sys.executable, str(generate_py), label, str(length)]
            cmd.extend(repr(freq) for freq in freqs)
            cmd.append("1")
            proc = _run(cmd)
            if proc.returncode != 0:
                sys.stderr.write(proc.stderr)
                _die(f"generate_random_dna_fasta.py failed for {label}")
            handle.write(proc.stdout)


def _frequencies(freqs_py: Path, fasta_path: Path, tsv_path: Path) -> list[tuple[str, int, list[float]]]:
    proc = _run([sys.executable, str(freqs_py), str(fasta_path)])
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        _die("dna_fasta_freqs.py failed")
    tsv_path.write_text(proc.stdout, encoding="utf-8", newline="\n")
    rows: list[tuple[str, int, list[float]]] = []
    for line_number, line in enumerate(proc.stdout.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) != 6:
            _die(f"freqs.tsv line {line_number} has {len(fields)} fields")
        label = fields[0]
        length = int(fields[1])
        reported = [float(field) for field in fields[2:]]
        rows.append((label, length, reported))
    return rows


def _check_rows(
    specs: list[tuple[str, int, list[float]]],
    rows: list[tuple[str, int, list[float]]],
    fasta_path: Path,
) -> None:
    if len(rows) != len(specs):
        _die(f"expected {len(specs)} frequency rows, got {len(rows)}")

    records: list[tuple[str, str]] = []

    def onseq(label: str, seq: str) -> None:
        records.append((label, seq))

    fasta.ReadSeqsOnSeq(str(fasta_path), onseq)
    if len(records) != len(specs):
        _die(f"expected {len(specs)} FASTA records, got {len(records)}")

    for spec, row, record in zip(specs, rows, records):
        label, length, freqs = spec
        row_label, row_length, reported = row
        rec_label, seq = record
        if row_label != label or rec_label != label:
            _die(f"label mismatch for {label}")
        if row_length != length or len(seq) != length:
            _die(f"length mismatch for {label}")
        total = sum(freqs)
        probs = [freq / total for freq in freqs]
        seq_u = seq.upper()
        tolerance = 1.0 / length + 0.0005
        for base, freq, prob, value in zip(_BASES, freqs, probs, reported):
            count = seq_u.count(base)
            if abs(count - prob * length) >= 1.0:
                _die(
                    f"{label} {base} count {count} is outside rounding error of {prob * length}"
                )
            if abs(value - freq) > tolerance + 1e-12:
                _die(
                    f"{label} {base} reported {value} vs requested {freq} (tolerance {tolerance})"
                )


def main() -> None:
    py_dir = Path(__file__).resolve().parent
    out_dir = py_dir.parent / "test_dna_freqs"
    out_dir.mkdir(parents=True, exist_ok=True)
    generate_py = py_dir / "generate_random_dna_fasta.py"
    freqs_py = py_dir / "dna_fasta_freqs.py"
    fasta_path = out_dir / "random_dna.fasta"
    tsv_path = out_dir / "freqs.tsv"

    _check_fatal(generate_py)
    specs = _draw_specs(random.Random(_SEED))
    _generate(generate_py, specs, fasta_path)
    rows = _frequencies(freqs_py, fasta_path, tsv_path)
    _check_rows(specs, rows, fasta_path)
    print(f"{_N} sequences: frequencies agree within rounding error")
    print("frequencies summing to 1.5 rejected")


if __name__ == "__main__":
    main()
