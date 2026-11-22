#!/usr/bin/env python3
"""Shuffle one sequence, fold it with RNAfold --circ, and report the MFE Z-score.

Mononucleotide shuffles keep the base composition. The original and each
permutation are folded as circular RNA. The report on stdout gives an ASCII
histogram, the rank of the original, and

    Z = (original MFE - mean) / sample standard deviation

of the permuted MFEs. Python's random generator is not gawk's rand(), so a
seed shared with bash/permute_fold.sh does not reproduce that script's
sequences.

    python3 py/permute_fold.py input.fa
        [--permutations N] [--seed S] [--threads T] [--bins B] [--barwidth W]
        [--permuted-fa PATH] [--fold PATH] [--mfe PATH]
"""

from __future__ import annotations

import argparse
import math
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
from typing import NoReturn

_MFE_AT_END = re.compile(r"\( *-?[0-9]+\.[0-9]+\)$")
_POSITIVE_INT = re.compile(r"[1-9][0-9]*")
_NONNEG_INT = re.compile(r"[0-9]+")


def _die(message: str) -> NoReturn:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def z_score(original_mfe: float, mean: float, stdev: float) -> float | None:
    """Return (original MFE - mean) / sample stdev, or None when stdev is not positive."""
    if not (
        math.isfinite(original_mfe)
        and math.isfinite(mean)
        and math.isfinite(stdev)
        and stdev > 0.0
    ):
        return None
    return (original_mfe - mean) / stdev


def read_first_fasta(path: str) -> tuple[str, str]:
    """Return (id, sequence) of the first FASTA record. Id is the first header token."""
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            seq_id, seq, nrec = _scan_fasta(handle)
    except (OSError, UnicodeError):
        _die(f"ERROR: cannot read '{path}'")
    if nrec < 1:
        _die(f"ERROR: no FASTA header found in '{path}'")
    if nrec != 1:
        print(
            f"WARNING: {nrec} records in '{path}'; using the first only",
            file=sys.stderr,
        )
    if seq == "":
        _die(f"ERROR: empty sequence in '{path}'")
    return seq_id, seq


def _scan_fasta(handle) -> tuple[str, str, int]:
    seq_id = ""
    parts: list[str] = []
    nrec = 0
    in_first = False
    for raw in handle:
        line = raw.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            nrec += 1
            if nrec == 1:
                header = line[1:].strip()
                seq_id = header.split()[0] if header else ""
                in_first = True
            else:
                in_first = False
            continue
        if in_first:
            parts.append(line.translate({ord(c): None for c in " \t"}))
    return seq_id, "".join(parts), nrec


def permute_sequences(seq: str, seq_id: str, n: int, seed: int) -> list[tuple[str, str]]:
    """Return n mononucleotide shuffles as (id, sequence), Fisher-Yates, composition preserved."""
    rng = random.Random(seed)
    width = len(str(n))
    records: list[tuple[str, str]] = []
    for perm in range(1, n + 1):
        chars = list(seq)
        for i in range(len(chars) - 1, 0, -1):
            j = rng.randrange(i + 1)
            chars[i], chars[j] = chars[j], chars[i]
        records.append((f"{seq_id}_perm_{perm:0{width}d}", "".join(chars)))
    return records


def format_fasta(records: list[tuple[str, str]]) -> str:
    return "".join(f">{label}\n{seq}\n" for label, seq in records)


def run_rnafold(fasta_text: str, threads: int) -> str:
    exe = shutil.which("RNAfold")
    if exe is None:
        _die("ERROR: RNAfold not found in PATH")
    try:
        proc = subprocess.run(
            [exe, "--circ", "--noPS", f"--jobs={threads}"],
            input=fasta_text,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
    except UnicodeError:
        _die("ERROR: RNAfold output was not valid UTF-8")
    if proc.returncode != 0:
        detail = proc.stderr.strip().splitlines()
        tail = detail[-1] if detail else ""
        message = f"ERROR: RNAfold exited {proc.returncode}"
        if tail:
            message = f"{message}: {tail}"
        _die(message)
    return proc.stdout


def parse_fold(text: str) -> tuple[str, list[tuple[str, str]]]:
    """Return (original structure, [(id, mfe text), ...]) in fold order.

    The structure is the first whitespace field of the original structure line.
    MFE text is the trailing parenthesized number, unchanged apart from
    removing parentheses and spaces.
    """
    seq_id = ""
    structure = ""
    rows: list[tuple[str, str]] = []
    for raw in text.splitlines():
        line = raw.rstrip("\r")
        if line == "":
            continue
        if line.startswith(">"):
            header = line[1:].strip()
            seq_id = header.split()[0] if header else ""
            continue
        found = _MFE_AT_END.search(line)
        if found is None:
            continue
        energy = re.sub(r"[() ]", "", found.group(0))
        if not rows:
            fields = line.split()
            structure = fields[0] if fields else ""
        rows.append((seq_id, energy))
    return structure, rows


def format_mfe_tsv(rows: list[tuple[str, str]]) -> str:
    return "".join(f"{label}\t{energy}\n" for label, energy in rows)


def _nice_bins(lo: float, hi: float, nbins: int) -> tuple[float, float, int]:
    """Return (start, width, n_bins) with the 1/2/5 x 10^k widths of permute_fold.sh."""
    raw = (hi - lo) / nbins
    if raw <= 0.0:
        raw = 0.1
    mag = 10.0 ** math.floor(math.log(raw) / math.log(10.0))
    if not math.isfinite(mag) or mag <= 0.0:
        width = 0.1
    else:
        fraction = raw / mag
        if fraction <= 1.0:
            step = 1.0
        elif fraction <= 2.0:
            step = 2.0
        elif fraction <= 5.0:
            step = 5.0
        else:
            step = 10.0
        width = step * mag
        if width < 0.1:
            width = 0.1
    start = math.floor(lo / width + 1e-9) * width
    n_bins = math.floor((hi - start) / width + 1e-9) + 1
    if n_bins < 1:
        n_bins = 1
    return start, width, n_bins


def _fold60(text: str) -> list[str]:
    return [f"    {text[i : i + 60]}" for i in range(0, len(text), 60)]


def format_report(
    *,
    infile: str,
    seq_id: str,
    seq: str,
    structure: str,
    seed: int,
    permuted_fa: str | None,
    original_mfe: float,
    permuted_mfe: list[float],
    nbins: int,
    barwidth: int,
) -> str:
    """ASCII histogram, rank, empirical p, and Z-score. Z is omitted when stdev is 0."""
    n = len(permuted_mfe)
    mean = statistics.mean(permuted_mfe)
    stdev = statistics.stdev(permuted_mfe) if n > 1 else 0.0
    min_mfe = min(permuted_mfe)
    max_mfe = max(permuted_mfe)
    lower = 0
    tied = 0
    for energy in permuted_mfe:
        if energy < original_mfe:
            lower += 1
        elif energy == original_mfe:
            tied += 1
    rank = lower + 1
    p_value = (lower + tied + 1) / (n + 1)
    z = z_score(original_mfe, mean, stdev)

    lo = original_mfe if original_mfe < min_mfe else min_mfe
    hi = original_mfe if original_mfe > max_mfe else max_mfe
    start, width, n_bins = _nice_bins(lo, hi, nbins)
    counts = [0] * n_bins
    for energy in permuted_mfe:
        index = math.floor((energy - start) / width + 1e-9)
        if 0 <= index < n_bins:
            counts[index] += 1
    original_bin = math.floor((original_mfe - start) / width + 1e-9)
    cmax = max(counts) if counts else 0
    scale = cmax / barwidth if cmax > barwidth else 1.0

    permfa_text = permuted_fa if permuted_fa else "(not written)"
    lines = [
        "=" * 78,
        f"RNAfold --circ MFE: original vs. {n} permuted sequences",
        "=" * 78,
        f"Input file     : {infile}",
        f"Sequence       : {seq_id} ({len(seq)} nt)",
        f"Permutations   : {n}, mononucleotide shuffle, seed {seed}",
        f"Permuted fasta : {permfa_text}",
        "",
        "Original sequence and circular MFE structure:",
    ]
    lines.extend(_fold60(seq))
    lines.append("")
    lines.extend(_fold60(structure))
    lines.append("")
    lines.append("Original MFE   : %8.2f kcal/mol" % original_mfe)
    lines.append(
        "Permuted MFE   : mean %.2f, sd %.2f, min %.2f, max %.2f kcal/mol"
        % (mean, stdev, min_mfe, max_mfe)
    )
    lines.append("")
    lines.append(
        "Histogram of permuted MFEs (bin width %g kcal/mol, each # ~ %.2g sequences)"
        % (width, scale)
    )
    lines.append("")
    lines.append("  %19s  %6s" % ("MFE bin (kcal/mol)", "count"))
    for bin_index, count in enumerate(counts):
        bar_len = int(count / scale + 0.5)
        if count > 0 and bar_len == 0:
            bar_len = 1
        mark = ""
        if bin_index == original_bin:
            mark = "  <== ORIGINAL (%.2f)" % original_mfe
        left = start + bin_index * width
        right = start + (bin_index + 1) * width
        lines.append(
            "  [%8.2f,%8.2f)  %6d  %s%s" % (left, right, count, "#" * bar_len, mark)
        )
    higher = n - lower - tied
    lines.append("")
    lines.append("Rank order (1 = lowest / most stable MFE, among original + permutations):")
    lines.append("  Original rank          : %d of %d" % (rank, n + 1))
    lines.append(
        "  Permuted with MFE <  original : %d (%.1f%%)" % (lower, 100.0 * lower / n)
    )
    lines.append("  Permuted with MFE == original : %d" % tied)
    lines.append(
        "  Permuted with MFE >  original : %d (%.1f%%)" % (higher, 100.0 * higher / n)
    )
    if z is not None:
        lines.append("  Z-score                : %.2f" % z)
    lines.append("  Empirical p (MFE <= original) : %.4g" % p_value)
    lines.append("=" * 78)
    return "\n".join(lines) + "\n"


def _write(path: str, text: str) -> None:
    try:
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
    except OSError as exc:
        detail = exc.strerror or "write failed"
        _die(f"ERROR: cannot write '{path}': {detail}")


def _positive_int(text: str, message: str) -> int:
    if _POSITIVE_INT.fullmatch(text) is None:
        _die(message)
    return int(text)


def _nonneg_int(text: str, message: str) -> int:
    if _NONNEG_INT.fullmatch(text) is None:
        _die(message)
    return int(text)


def _default_threads() -> int:
    count = os.cpu_count() or 1
    return max(1, min(count, 20))


def _from_env_or_flag(flag: str | None, env_name: str, fallback: int, message: str) -> int:
    if flag is not None:
        return _positive_int(flag, message)
    env = os.environ.get(env_name)
    if env is not None and env.strip() != "":
        return _positive_int(env.strip(), message)
    return fallback


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Shuffle a sequence, fold the original and each permutation with "
            "RNAfold --circ, and print an MFE histogram, rank, and Z-score."
        )
    )
    parser.add_argument("input_fa", metavar="input.fa", help="Input FASTA; the first record is folded")
    parser.add_argument(
        "--permutations",
        default="1000",
        help="Number of mononucleotide shuffles (default: 1000)",
    )
    parser.add_argument("--seed", default="1", help="Non-negative integer seed (default: 1)")
    parser.add_argument(
        "--threads",
        default=None,
        help="RNAfold --jobs (default: $THREADS, else min(cpu count, 20))",
    )
    parser.add_argument(
        "--bins",
        default=None,
        help="Approximate histogram bins (default: $BINS, else 20)",
    )
    parser.add_argument(
        "--barwidth",
        default=None,
        help="Width of the longest histogram bar (default: $BARWIDTH, else 60)",
    )
    parser.add_argument(
        "--permuted-fa",
        default=None,
        metavar="PATH",
        help="Write shuffled sequences to PATH (default: do not write)",
    )
    parser.add_argument(
        "--fold",
        default=None,
        metavar="PATH",
        help="Write raw RNAfold output to PATH (default: do not write)",
    )
    parser.add_argument(
        "--mfe",
        default=None,
        metavar="PATH",
        help="Write id and MFE TSV to PATH (default: do not write)",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    infile = args.input_fa
    # Match permute_fold.sh: a missing file is reported before the other checks.
    if not os.path.isfile(infile) or not os.access(infile, os.R_OK):
        _die(f"ERROR: cannot read '{infile}'")
    n = _positive_int(args.permutations, "ERROR: n_permutations must be a positive integer")
    seed = _nonneg_int(args.seed, "ERROR: seed must be a non-negative integer")
    threads = _from_env_or_flag(
        args.threads,
        "THREADS",
        _default_threads(),
        "ERROR: threads must be a positive integer",
    )
    nbins = _from_env_or_flag(args.bins, "BINS", 20, "ERROR: bins must be a positive integer")
    barwidth = _from_env_or_flag(
        args.barwidth,
        "BARWIDTH",
        60,
        "ERROR: barwidth must be a positive integer",
    )
    if shutil.which("RNAfold") is None:
        _die("ERROR: RNAfold not found in PATH")

    seq_id, seq = read_first_fasta(infile)
    permuted = permute_sequences(seq, seq_id, n, seed)
    permuted_text = format_fasta(permuted)
    if args.permuted_fa:
        _write(args.permuted_fa, permuted_text)

    fold_text = run_rnafold(format_fasta([(seq_id, seq)]) + permuted_text, threads)
    if args.fold:
        _write(args.fold, fold_text)

    structure, rows = parse_fold(fold_text)
    if args.mfe:
        _write(args.mfe, format_mfe_tsv(rows))
    if len(rows) != n + 1:
        source = args.fold if args.fold else "RNAfold output"
        _die(f"ERROR: expected {n + 1} MFE values, parsed {len(rows)} from {source}")

    original_mfe = float(rows[0][1])
    permuted_mfe = [float(energy) for _, energy in rows[1:]]
    report = format_report(
        infile=infile,
        seq_id=seq_id,
        seq=seq,
        structure=structure,
        seed=seed,
        permuted_fa=args.permuted_fa,
        original_mfe=original_mfe,
        permuted_mfe=permuted_mfe,
        nbins=nbins,
        barwidth=barwidth,
    )
    sys.stdout.write(report)


if __name__ == "__main__":
    main()
