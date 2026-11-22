#!/usr/bin/env python3
"""Correlate measured and predicted circular-RNA MFE Z-scores.

Measured Z follows permute_fold: mononucleotide shuffles, RNAfold --circ, and

    Z = (original MFE - mean) / sample standard deviation

Predicted mean and standard deviation come from mfe_grid_interp. The grid
stores cmean and cstddev as 100 * mean / L and 100 * stddev / L, so they are
scaled back by L/100 before the same Z formula is applied.

    python3 py/correlate_measured_predicted_Z.py
        [--input DIR|FASTA] [--table PATH] [--outdir DIR]
        [--permutations N] [--threads T] [--seed S] [--limit N] [--sample N]
        [--tsv NAME] [--compare CONTROLS.tsv]
"""

from __future__ import annotations

import argparse
import io
import math
import os
import shutil
import statistics
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import mfe_grid_interp
import permute_fold

_ROOT = Path(__file__).resolve().parent.parent
_COLUMNS = (
    "label",
    "L",
    "fA",
    "fC",
    "fG",
    "fT",
    "MFE",
    "measured_mean",
    "measured_sd",
    "measured_Z",
    "pred_mean",
    "pred_sd",
    "predicted_Z",
    "status",
)
_FLOAT_COLUMNS = (
    "fA",
    "fC",
    "fG",
    "fT",
    "MFE",
    "measured_mean",
    "measured_sd",
    "measured_Z",
    "pred_mean",
    "pred_sd",
    "predicted_Z",
)


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _blank_row(label: str) -> dict[str, object]:
    row: dict[str, object] = {name: None for name in _COLUMNS}
    row["label"] = label
    row["status"] = ""
    return row


def _fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.6g}"
    return str(value)


def _format_row(row: dict[str, object]) -> str:
    return "\t".join(_fmt(row[name]) for name in _COLUMNS)


def _capture_system_exit(fn):
    """Run fn. On SystemExit, return (None, message); otherwise (value, None)."""
    captured = io.StringIO()
    old = sys.stderr
    sys.stderr = captured
    try:
        try:
            return fn(), None
        except SystemExit:
            message = " ".join(captured.getvalue().split()) or "failed"
            return None, message
    finally:
        sys.stderr = old


def _base_counts(seq: str) -> np.ndarray | None:
    """Return ACGT counts, or None when seq contains any other base."""
    counts = np.zeros(4, dtype=np.float64)
    index = {"A": 0, "C": 1, "G": 2, "T": 3}
    for base in seq:
        slot = index.get(base)
        if slot is None:
            return None
        counts[slot] += 1.0
    return counts


def _in_length_range(lengths: np.ndarray, length: int) -> bool:
    return float(lengths[0]) - 1e-6 <= length <= float(lengths[-1]) + 1e-6


def _prepare_sequence(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _sequence_can_score(seq: str, lengths: np.ndarray) -> bool:
    """True when the sequence is ACGT and inside the tabulated length range."""
    seq = _prepare_sequence(seq)
    if not seq or any(base not in "ACGT" for base in seq):
        return False
    return _in_length_range(lengths, len(seq))


def _can_score(path: Path, lengths: np.ndarray) -> bool:
    """True when the file's first sequence is ACGT and inside the tabulated length range."""
    loaded, error = _capture_system_exit(lambda: permute_fold.read_first_fasta(str(path)))
    if error is not None or loaded is None:
        return False
    _seq_id, seq = loaded
    return _sequence_can_score(seq, lengths)


def _evenly_spaced(items: list, count: int) -> list:
    """Take count items spread across a list, including both ends."""
    if count >= len(items):
        return list(items)
    last = len(items) - 1
    return [items[int(round(i * last / (count - 1)))] for i in range(count)]


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
    except (OSError, UnicodeError):
        _die(f"ERROR: cannot read '{path}'")
    if seq_id or parts:
        finish()
    return records


def score_sequence(
    seq_id: str,
    seq: str,
    fallback_label: str,
    table_rows: dict,
    lengths: np.ndarray,
    k: int,
    n_permutations: int,
    seed: int,
    threads: int,
) -> tuple[dict[str, object], str]:
    """Return (row, detail). detail is extra text for the progress line."""
    label = seq_id or fallback_label
    label = label.replace("\t", " ").replace("\r", " ").replace("\n", " ")
    row = _blank_row(label)
    seq = _prepare_sequence(seq)
    length = len(seq)
    row["L"] = length
    counts = _base_counts(seq)
    if counts is None or length <= 0:
        row["status"] = "non_acgt"
        return row, ""
    freqs = counts / float(length)
    row["fA"] = float(freqs[0])
    row["fC"] = float(freqs[1])
    row["fG"] = float(freqs[2])
    row["fT"] = float(freqs[3])
    if not _in_length_range(lengths, length):
        row["status"] = "length_out_of_range"
        return row, ""

    interpolated, error = _capture_system_exit(
        lambda: mfe_grid_interp.interpolate(table_rows, lengths, k, float(length), freqs)
    )
    if error is not None or interpolated is None:
        row["status"] = "interp_failed"
        return row, error or "interp_failed"
    cmean, cstddev, _support = interpolated
    pred_mean = float(cmean) * length / 100.0
    pred_sd = float(cstddev) * length / 100.0
    row["pred_mean"] = pred_mean
    row["pred_sd"] = pred_sd
    if not (math.isfinite(pred_sd) and pred_sd > 0.0):
        row["status"] = "zero_pred_sd"
        return row, ""

    def fold():
        permuted = permute_fold.permute_sequences(seq, seq_id, n_permutations, seed)
        fasta_text = permute_fold.format_fasta([(seq_id, seq)])
        fasta_text += permute_fold.format_fasta(permuted)
        fold_text = permute_fold.run_rnafold(fasta_text, threads)
        return permute_fold.parse_fold(fold_text)

    parsed, error = _capture_system_exit(fold)
    if error is not None or parsed is None:
        row["status"] = "fold_failed"
        return row, error or "fold_failed"
    _structure, mfe_rows = parsed
    if len(mfe_rows) != n_permutations + 1:
        row["status"] = "fold_failed"
        return row, f"expected {n_permutations + 1} MFE values, parsed {len(mfe_rows)}"
    try:
        original_mfe = float(mfe_rows[0][1])
        permuted_mfe = [float(energy) for _label, energy in mfe_rows[1:]]
    except ValueError as exc:
        row["status"] = "fold_failed"
        return row, str(exc)
    measured_mean = statistics.mean(permuted_mfe)
    measured_sd = statistics.stdev(permuted_mfe) if len(permuted_mfe) > 1 else 0.0
    row["MFE"] = original_mfe
    row["measured_mean"] = measured_mean
    row["measured_sd"] = measured_sd
    row["predicted_Z"] = permute_fold.z_score(original_mfe, pred_mean, pred_sd)
    row["measured_Z"] = permute_fold.z_score(original_mfe, measured_mean, measured_sd)
    if row["measured_Z"] is None:
        row["status"] = "zero_measured_sd"
        return row, ""
    row["status"] = "ok"
    return row, ""


def score_file(
    path: Path,
    table_rows: dict,
    lengths: np.ndarray,
    k: int,
    n_permutations: int,
    seed: int,
    threads: int,
) -> tuple[dict[str, object], str]:
    """Score the first FASTA record in path."""
    loaded, error = _capture_system_exit(lambda: permute_fold.read_first_fasta(str(path)))
    if error is not None or loaded is None:
        row = _blank_row(path.name)
        row["status"] = "read_failed"
        return row, error or "read_failed"
    seq_id, seq = loaded
    return score_sequence(
        seq_id,
        seq,
        path.name,
        table_rows,
        lengths,
        k,
        n_permutations,
        seed,
        threads,
    )


def _pearson(measured: list[float], predicted: list[float]) -> float:
    if len(measured) < 2:
        return float("nan")
    coefficient = np.corrcoef(np.asarray(measured), np.asarray(predicted))[0, 1]
    return float(coefficient)


def _write_scatter(path: Path, measured: list[float], predicted: list[float], r: float) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(measured, predicted, s=12, alpha=0.75)
    ax.set_xlabel("Measured Z")
    ax.set_ylabel("Predicted Z")
    ax.text(
        0.04,
        0.96,
        "r = %.3f\nn = %d" % (r, len(measured)),
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _pairs_from_tsv(path: Path) -> tuple[list[float], list[float], list[str]]:
    """Return finite (measured Z, predicted Z, label) rows from a score table."""
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            try:
                label_col = header.index("label")
                measured_col = header.index("measured_Z")
                predicted_col = header.index("predicted_Z")
            except ValueError:
                _die(f"ERROR: {path} is missing label, measured_Z, or predicted_Z")
            labels: list[str] = []
            measured: list[float] = []
            predicted: list[float] = []
            width = max(label_col, measured_col, predicted_col)
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) <= width:
                    continue
                try:
                    measured_z = float(fields[measured_col])
                    predicted_z = float(fields[predicted_col])
                except ValueError:
                    continue
                if math.isfinite(measured_z) and math.isfinite(predicted_z):
                    labels.append(fields[label_col])
                    measured.append(measured_z)
                    predicted.append(predicted_z)
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    return measured, predicted, labels


def _control_points(
    measured: list[float], predicted: list[float], labels: list[str], prefix: str
) -> tuple[list[float], list[float]]:
    xs = []
    ys = []
    for x, y, label in zip(measured, predicted, labels):
        if label.startswith(prefix):
            xs.append(x)
            ys.append(y)
    return xs, ys


def _write_combined(
    path: Path,
    controls: tuple[list[float], list[float], list[str]],
    rods: tuple[list[float], list[float], list[str]],
    limits: tuple[float, float] | None = None,
) -> tuple[float, float, float]:
    """Scatter DECOY, VLE, and RODs. Return Pearson r for RODs, controls, and both.

    When limits is set, both axes use that range so the dashed y=x line is the
    diagonal of the axes box.
    """
    control_measured, control_predicted, control_labels = controls
    rod_measured, rod_predicted, _rod_labels = rods
    r_rods = _pearson(rod_measured, rod_predicted)
    r_controls = _pearson(control_measured, control_predicted)
    r_both = _pearson(control_measured + rod_measured, control_predicted + rod_predicted)
    decoy_x, decoy_y = _control_points(control_measured, control_predicted, control_labels, "DECOY")
    vle_x, vle_y = _control_points(control_measured, control_predicted, control_labels, "VLE")
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(decoy_x, decoy_y, s=12, alpha=0.75, color="orange", label="DECOY")
    ax.scatter(vle_x, vle_y, s=28, alpha=0.9, color="magenta", label="VLE")
    ax.scatter(rod_measured, rod_predicted, s=12, alpha=0.75, color="cyan", label="RODs")
    ax.axline((0, 0), slope=1, linestyle="--", color="0.35", linewidth=1, zorder=0)
    if limits is not None:
        low, high = limits
        ax.set_xlim(low, high)
        ax.set_ylim(low, high)
        ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Measured Z")
    ax.set_ylabel("Predicted Z")
    ax.legend(loc="lower right")
    ax.text(
        0.04,
        0.96,
        "RODs r = %.3f (n = %d)\ncontrols r = %.3f (n = %d)\nboth r = %.3f (n = %d)"
        % (
            r_rods,
            len(rod_measured),
            r_controls,
            len(control_measured),
            r_both,
            len(control_measured) + len(rod_measured),
        ),
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return r_rods, r_controls, r_both


def _progress(index: int, total: int, row: dict[str, object], detail: str) -> None:
    fields = [
        f"{index}/{total}",
        str(row["label"]),
        str(row["status"]),
        "L=%s" % _fmt(row["L"]),
        "measured_Z=%s" % _fmt(row["measured_Z"]),
        "predicted_Z=%s" % _fmt(row["predicted_Z"]),
    ]
    if detail:
        fields.append(detail)
    print("\t".join(fields), file=sys.stderr, flush=True)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Correlate permutation MFE Z-scores with Z-scores predicted by "
            "interpolating the RNAfold length and composition grid."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=_ROOT / "controls1000",
        help="Directory of one-sequence FASTA files, or one multi-FASTA file",
    )
    parser.add_argument(
        "--table",
        type=Path,
        default=_ROOT / "output" / "rnafold_grid_mean_var.tsv",
        help="RNAfold grid mean/variance table",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=_ROOT / "correlate_measured_predicted_Z",
        help="Directory for the TSV and scatterplot",
    )
    parser.add_argument("--permutations", type=int, default=100, help="Mononucleotide shuffles per sequence")
    parser.add_argument("--threads", type=int, default=64, help="Maximum RNAfold --jobs")
    parser.add_argument("--seed", type=int, default=1, help="Shuffle seed")
    parser.add_argument("--limit", type=int, default=0, help="Score only the first N files (0 = all)")
    parser.add_argument(
        "--sample",
        type=int,
        default=0,
        help="Evenly sample N in-range ACGT sequences from the sorted input (0 = no sampling)",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=Path("z_scores.tsv"),
        help="Score table name, relative to --outdir unless absolute",
    )
    parser.add_argument(
        "--compare",
        type=Path,
        default=None,
        help="Control score table to plot with this run; writes scatter.png",
    )
    args = parser.parse_args(argv)

    if args.permutations < 1:
        _die("ERROR: permutations must be a positive integer")
    if args.threads < 1:
        _die("ERROR: threads must be a positive integer")
    if args.seed < 0:
        _die("ERROR: seed must be a non-negative integer")
    if args.limit < 0:
        _die("ERROR: limit must be a non-negative integer")
    if args.sample < 0:
        _die("ERROR: sample must be a non-negative integer")
    if args.sample and args.sample < 2:
        _die("ERROR: sample must be 0 or at least 2")
    if not args.input.exists():
        _die(f"ERROR: not found: {args.input}")
    if shutil.which("RNAfold") is None:
        _die("ERROR: RNAfold not found in PATH")

    table_rows, lengths, k = mfe_grid_interp.load_table(str(args.table))
    records: list[tuple[str, str]] | None = None
    files: list[Path] = []
    if args.input.is_dir():
        files = sorted(path for path in args.input.iterdir() if path.is_file())
        n_input = len(files)
        if args.sample:
            candidates = [path for path in files if _can_score(path, lengths)]
            if len(candidates) < 2:
                _die(f"ERROR: only {len(candidates)} scorable sequences in {args.input}")
            files = _evenly_spaced(candidates, args.sample)
            print(
                f"sample {len(files)} of {len(candidates)} scorable from {n_input} files",
                file=sys.stderr,
                flush=True,
            )
        elif args.limit:
            files = files[: args.limit]
        if not files:
            _die(f"ERROR: no files in {args.input}")
    elif args.input.is_file():
        records = _read_fasta(args.input)
        n_input = len(records)
        if args.sample:
            candidates = [record for record in records if _sequence_can_score(record[1], lengths)]
            if len(candidates) < 2:
                _die(f"ERROR: only {len(candidates)} scorable sequences in {args.input}")
            records = _evenly_spaced(candidates, args.sample)
            print(
                f"sample {len(records)} of {len(candidates)} scorable from {n_input} records",
                file=sys.stderr,
                flush=True,
            )
        elif args.limit:
            records = records[: args.limit]
        if not records:
            _die(f"ERROR: no FASTA records in {args.input}")
    else:
        _die(f"ERROR: not a file or directory: {args.input}")

    cpu = os.cpu_count() or 1
    threads = max(1, min(args.threads, cpu, args.permutations + 1))
    args.outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = args.tsv if args.tsv.is_absolute() else args.outdir / args.tsv
    png_path = args.outdir / "scatter.png"
    compare_path = None
    if args.compare is not None:
        compare_path = args.compare if args.compare.is_absolute() else args.outdir / args.compare
        if not compare_path.is_file():
            _die(f"ERROR: comparison table not found: {compare_path}")
    print(
        f"sequences {len(records) if records is not None else len(files)} "
        f"permutations {args.permutations} RNAfold_jobs {threads}",
        file=sys.stderr,
        flush=True,
    )

    measured: list[float] = []
    predicted: list[float] = []
    with tsv_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(_COLUMNS) + "\n")
        if records is not None:
            jobs = [
                (index, None, seq_id, seq)
                for index, (seq_id, seq) in enumerate(records, start=1)
            ]
        else:
            jobs = [(index, path, "", "") for index, path in enumerate(files, start=1)]
        for index, path, seq_id, seq in jobs:
            if path is None:
                row, detail = score_sequence(
                    seq_id,
                    seq,
                    seq_id or f"record_{index}",
                    table_rows,
                    lengths,
                    k,
                    args.permutations,
                    args.seed,
                    threads,
                )
            else:
                row, detail = score_file(
                    path,
                    table_rows,
                    lengths,
                    k,
                    args.permutations,
                    args.seed,
                    threads,
                )
            handle.write(_format_row(row) + "\n")
            handle.flush()
            total = len(records) if records is not None else len(files)
            _progress(index, total, row, detail)
            measured_z = row["measured_Z"]
            predicted_z = row["predicted_Z"]
            if (
                isinstance(measured_z, float)
                and isinstance(predicted_z, float)
                and math.isfinite(measured_z)
                and math.isfinite(predicted_z)
            ):
                measured.append(measured_z)
                predicted.append(predicted_z)

    if compare_path is not None:
        rods = _pairs_from_tsv(tsv_path)
        controls = _pairs_from_tsv(compare_path)
        r_rods, r_controls, r_both = _write_combined(png_path, controls, rods)
        print("Pearson r RODs = %.6g (n = %d)" % (r_rods, len(rods[0])))
        print("Pearson r controls = %.6g (n = %d)" % (r_controls, len(controls[0])))
        print("Pearson r both = %.6g (n = %d)" % (r_both, len(rods[0]) + len(controls[0])))
    else:
        r = _pearson(measured, predicted)
        if len(measured) >= 2:
            _write_scatter(png_path, measured, predicted, r)
        elif png_path.exists():
            png_path.unlink()
        print("Pearson r = %.6g (n = %d)" % (r, len(measured)))


if __name__ == "__main__":
    main()
