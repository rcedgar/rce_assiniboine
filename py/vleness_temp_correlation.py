#!/usr/bin/env python3
"""Pearson correlations of vleness S across folding temperatures.

Rows are the labels in joined_T0_T10_T37_T60.tsv. S is read from
T0/T10/T37/T60.vleness.tsv. Writes report.txt, a six-panel scatter of S,
histograms of S(T) - S(T37), octave histograms of dot-bracket Hamming
distance versus T37, and a four-panel histogram of S(T0) - S(T37) split
into VLE and DECOY into big_rnafold_controls_100_per_family_temps/.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
_DIR = _ROOT / "big_rnafold_controls_100_per_family_temps"
_TEMPS = ("T0", "T10", "T37", "T60")
_PAIRS = (
    ("T0", "T10"),
    ("T0", "T37"),
    ("T0", "T60"),
    ("T10", "T37"),
    ("T10", "T60"),
    ("T37", "T60"),
)
_VERSUS_T37 = ("T0", "T10", "T60")
_SS_COLUMN = {"T0": 2, "T10": 4, "T37": 6, "T60": 8}
_N_DELTA_BINS = 33
_JOINED_COLUMNS = 9
_CLASSES = ("VLE", "DECOY")
_CLASS_COLOR = {"VLE": "#d55e00", "DECOY": "#0072b2"}


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _save_fig(fig: plt.Figure, path: Path) -> bool:
    try:
        fig.savefig(path, dpi=150)
    except OSError as exc:
        print("ERROR: cannot write '%s': %s" % (path, exc.strerror or exc), file=sys.stderr)
        plt.close(fig)
        return False
    plt.close(fig)
    print("wrote %s" % path)
    return True


def _load_scores(path: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            try:
                label_col = header.index("label")
                score_col = header.index("S")
            except ValueError:
                _die("ERROR: %s is missing label or S" % path)
            width = max(label_col, score_col)
            for line_number, line in enumerate(handle, start=2):
                text = line.rstrip("\r\n")
                if not text:
                    continue
                fields = text.split("\t")
                if len(fields) <= width or fields[score_col] == "":
                    _die("ERROR: %s:%d is missing S" % (path, line_number))
                label = fields[label_col]
                if label in scores:
                    _die("ERROR: %s:%d repeats label %s" % (path, line_number, label))
                try:
                    value = float(fields[score_col])
                except ValueError:
                    _die("ERROR: %s:%d S is not a number" % (path, line_number))
                if not math.isfinite(value):
                    _die("ERROR: %s:%d S is not finite" % (path, line_number))
                scores[label] = value
    except OSError as exc:
        _die("ERROR: cannot read '%s': %s" % (path, exc.strerror or exc))
    if not scores:
        _die("ERROR: no S values in %s" % path)
    return scores


def _load_joined(path: Path) -> list[tuple[str, dict[str, str]]]:
    rows: list[tuple[str, dict[str, str]]] = []
    seen: set[str] = set()
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for line_number, line in enumerate(handle, start=1):
                text = line.rstrip("\r\n")
                if not text:
                    continue
                fields = text.split("\t")
                if len(fields) != _JOINED_COLUMNS:
                    _die(
                        "ERROR: %s:%d has %d columns"
                        % (path, line_number, len(fields))
                    )
                label = fields[0]
                if label in seen:
                    _die("ERROR: %s:%d repeats label %s" % (path, line_number, label))
                seen.add(label)
                structures = {temp: fields[column] for temp, column in _SS_COLUMN.items()}
                baseline = structures["T37"]
                for temp, structure in structures.items():
                    if len(structure) != len(baseline):
                        _die(
                            "ERROR: %s:%d %s length %d differs from T37 length %d"
                            % (path, line_number, temp, len(structure), len(baseline))
                        )
                rows.append((label, structures))
    except OSError as exc:
        _die("ERROR: cannot read '%s': %s" % (path, exc.strerror or exc))
    if not rows:
        _die("ERROR: no rows in %s" % path)
    return rows


def _scores_for_labels(labels: list[str], path: Path) -> np.ndarray:
    table = _load_scores(path)
    missing = [label for label in labels if label not in table]
    if missing:
        _die(
            "ERROR: %d joined labels missing S in %s, first %s"
            % (len(missing), path.name, missing[0])
        )
    return np.asarray([table[label] for label in labels], dtype=np.float64)


def _class_name(label: str) -> str:
    if label.startswith("VLE_"):
        return "VLE"
    if label.startswith("DECOY_"):
        return "DECOY"
    _die("ERROR: label is not VLE_ or DECOY_: %s" % label)


def _hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))


def _octave_index(count: int) -> int:
    if count == 0:
        return 0
    if count < 0:
        _die("ERROR: negative dot-bracket difference")
    return count.bit_length()


def _octave_bins(max_count: int) -> list[tuple[int, int, str]]:
    """Bins are 0, then 1, 2-3, 4-7, through the octave that holds max_count."""
    bins: list[tuple[int, int, str]] = [(0, 0, "0")]
    if max_count > 0:
        k = 1
        while True:
            lo = 1 << (k - 1)
            hi = (1 << k) - 1
            label = str(lo) if lo == hi else "%d-%d" % (lo, hi)
            bins.append((lo, hi, label))
            if hi >= max_count:
                break
            k += 1
    for index, (lo, hi, label) in enumerate(bins):
        if _octave_index(lo) != index or _octave_index(hi) != index:
            _die("ERROR: octave bin %s is misaligned" % label)
    return bins


def _bin_counts(values: np.ndarray, n_bins: int) -> np.ndarray:
    counts = np.zeros(n_bins, dtype=np.int64)
    for value in values:
        index = _octave_index(int(value))
        if index >= n_bins:
            _die("ERROR: difference %d is outside the octave bins" % int(value))
        counts[index] += 1
    if int(counts.sum()) != values.size:
        _die("ERROR: octave bins do not cover every row")
    return counts


def _mean_sd(values: np.ndarray) -> tuple[float, float]:
    if values.size < 2:
        _die("ERROR: need at least 2 values for a standard deviation")
    return float(values.mean()), float(values.std(ddof=1))


def _write_scatter(
    path: Path,
    scores: dict[str, np.ndarray],
    correlations: list[tuple[str, str, float]],
) -> bool:
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for ax, (left, right, pearson) in zip(axes.ravel(), correlations):
        x = scores[left]
        y = scores[right]
        low = float(min(x.min(), y.min()))
        high = float(max(x.max(), y.max()))
        pad = 0.05 * (high - low) if high > low else 1.0
        low -= pad
        high += pad
        ax.scatter(x, y, s=10, alpha=0.4, color="cyan", linewidths=0, zorder=2)
        ax.axline((0, 0), slope=1, linestyle="--", color="0.35", linewidth=1, zorder=1)
        ax.set_xlim(low, high)
        ax.set_ylim(low, high)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("S (%s)" % left)
        ax.set_ylabel("S (%s)" % right)
        ax.set_title("%s vs %s  r = %.3f" % (left, right, pearson))
    fig.tight_layout()
    return _save_fig(fig, path)


def _write_delta_histograms(path: Path, deltas: dict[str, np.ndarray]) -> bool:
    combined = np.concatenate([deltas[temp] for temp in _VERSUS_T37])
    low = float(combined.min())
    high = float(combined.max())
    if high <= low:
        high = low + 1.0
    edges = np.linspace(low, high, _N_DELTA_BINS + 1)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True, sharey=True)
    for ax, temp in zip(axes, _VERSUS_T37):
        ax.hist(deltas[temp], bins=edges, color="#4c78a8")
        ax.axvline(0.0, color="0.25", linestyle="--", linewidth=1)
        ax.set_title("%s - T37" % temp)
        ax.set_xlabel("S(%s) - S(T37)" % temp)
        ax.set_ylabel("Count")
    fig.tight_layout()
    return _save_fig(fig, path)


def _write_class_delta_histograms(path: Path, delta: np.ndarray, classes: np.ndarray) -> bool:
    """Counts on the top row, fraction of each class on the bottom. Shared bins."""
    series = {}
    for name in _CLASSES:
        values = delta[classes == name]
        if values.size == 0:
            _die("ERROR: no %s rows for S(T0) - S(T37)" % name)
        series[name] = values
    low = float(delta.min())
    high = float(delta.max())
    if high <= low:
        high = low + 1.0
    edges = np.linspace(low, high, _N_DELTA_BINS + 1)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    axes[0, 1].sharey(axes[0, 0])
    axes[1, 1].sharey(axes[1, 0])
    for column, name in enumerate(_CLASSES):
        values = series[name]
        color = _CLASS_COLOR[name]
        count_axis = axes[0, column]
        count_axis.hist(values, bins=edges, color=color)
        count_axis.axvline(0.0, color="0.25", linestyle="--", linewidth=1)
        count_axis.set_title("%s count (n = %d)" % (name, values.size))
        count_axis.set_ylabel("Count")
        fraction_axis = axes[1, column]
        weights = np.full(values.size, 1.0 / float(values.size))
        fraction_axis.hist(values, bins=edges, weights=weights, color=color)
        fraction_axis.axvline(0.0, color="0.25", linestyle="--", linewidth=1)
        fraction_axis.set_title("%s fraction (n = %d)" % (name, values.size))
        fraction_axis.set_xlabel("S(T0) - S(T37)")
        fraction_axis.set_ylabel("Fraction")
    fig.tight_layout()
    return _save_fig(fig, path)


def _write_diff_histograms(
    path: Path,
    bin_counts: dict[str, np.ndarray],
    bins: list[tuple[int, int, str]],
) -> bool:
    labels = [label for _lo, _hi, label in bins]
    positions = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True)
    for ax, temp in zip(axes, _VERSUS_T37):
        ax.bar(positions, bin_counts[temp], width=0.9, color="#4c78a8")
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_title("%s vs T37" % temp)
        ax.set_xlabel("Positions differing from T37")
        ax.set_ylabel("Count")
    fig.tight_layout()
    return _save_fig(fig, path)


def _write_report(
    path: Path,
    n: int,
    correlations: list[tuple[str, str, float]],
    scores: dict[str, np.ndarray],
    deltas: dict[str, np.ndarray],
    diffs: dict[str, np.ndarray],
    bin_counts: dict[str, np.ndarray],
    bins: list[tuple[int, int, str]],
    t60_matches_t37: bool,
    class_deltas: dict[str, np.ndarray],
) -> None:
    lines = ["n\t%d" % n, ""]
    if t60_matches_t37:
        lines.append(
            "T60 dot-brackets match T37 for every joined label, so Pearson r "
            "for T37 and T60 is 1 and every T60 minus T37 vleness change and "
            "Hamming distance is 0."
        )
        lines.append("")
    lines.append("pearson_r")
    lines.append("temp_a\ttemp_b\tr")
    for left, right, pearson in correlations:
        lines.append("%s\t%s\t%.6f" % (left, right, pearson))
    lines.append("")
    lines.append("S")
    lines.append("temp\tmean\tsd")
    for temp in _TEMPS:
        mean, sd = _mean_sd(scores[temp])
        lines.append("%s\t%.6f\t%.6f" % (temp, mean, sd))
    lines.append("")
    lines.append("delta_S_vs_T37")
    lines.append("temp\tmean\tsd")
    for temp in _VERSUS_T37:
        mean, sd = _mean_sd(deltas[temp])
        lines.append("%s\t%.6f\t%.6f" % (temp, mean, sd))
    lines.append("")
    lines.append("delta_S_T0_minus_T37_by_class")
    lines.append("class\tn\tmean\tsd")
    for name in _CLASSES:
        values = class_deltas[name]
        mean, sd = _mean_sd(values)
        lines.append("%s\t%d\t%.6f\t%.6f" % (name, values.size, mean, sd))
    lines.append("")
    lines.append("dotbracket_hamming_vs_T37")
    lines.append("temp\tn_identical\tmean\tmedian\tmax")
    for temp in _VERSUS_T37:
        values = diffs[temp]
        lines.append(
            "%s\t%d\t%.4f\t%.4f\t%d"
            % (
                temp,
                int(np.count_nonzero(values == 0)),
                float(values.mean()),
                float(np.median(values)),
                int(values.max()),
            )
        )
    lines.append("")
    lines.append("dotbracket_octave_bins")
    lines.append("bin\t" + "\t".join(_VERSUS_T37))
    for index, (_lo, _hi, label) in enumerate(bins):
        counts = "\t".join(str(int(bin_counts[temp][index])) for temp in _VERSUS_T37)
        lines.append("%s\t%s" % (label, counts))
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def main() -> None:
    joined_path = _DIR / "joined_T0_T10_T37_T60.tsv"
    rows = _load_joined(joined_path)
    labels = [label for label, _structures in rows]
    scores = {
        temp: _scores_for_labels(labels, _DIR / ("%s.vleness.tsv" % temp))
        for temp in _TEMPS
    }
    correlations = [
        (left, right, float(np.corrcoef(scores[left], scores[right])[0, 1]))
        for left, right in _PAIRS
    ]
    deltas = {temp: scores[temp] - scores["T37"] for temp in _VERSUS_T37}
    classes = np.asarray([_class_name(label) for label in labels])
    class_deltas = {name: deltas["T0"][classes == name] for name in _CLASSES}
    diffs = {
        temp: np.asarray(
            [
                _hamming(structures[temp], structures["T37"])
                for _label, structures in rows
            ],
            dtype=np.int64,
        )
        for temp in _VERSUS_T37
    }
    max_count = max(int(diffs[temp].max()) for temp in _VERSUS_T37)
    bins = _octave_bins(max_count)
    counts = {temp: _bin_counts(diffs[temp], len(bins)) for temp in _VERSUS_T37}
    t60_matches_t37 = bool(
        np.array_equal(scores["T37"], scores["T60"]) and np.all(diffs["T60"] == 0)
    )

    report_path = _DIR / "report.txt"
    scatter_path = _DIR / "vleness_pairs_scatter.png"
    delta_path = _DIR / "vleness_delta_T37_hist.png"
    diff_path = _DIR / "dotbracket_diffs_T37_hist.png"
    _write_report(
        report_path,
        len(labels),
        correlations,
        scores,
        deltas,
        diffs,
        counts,
        bins,
        t60_matches_t37,
        class_deltas,
    )
    class_delta_path = _DIR / "vleness_delta_T0_T37_by_class.png"
    wrote = [
        _write_class_delta_histograms(class_delta_path, deltas["T0"], classes),
        _write_scatter(scatter_path, scores, correlations),
        _write_delta_histograms(delta_path, deltas),
        _write_diff_histograms(diff_path, counts, bins),
    ]
    print("n = %d" % len(labels))
    for left, right, pearson in correlations:
        print("r %s %s = %.6f" % (left, right, pearson))
    print("wrote %s" % report_path)
    if not all(wrote):
        _die("ERROR: one or more figures could not be written")
    for name in _CLASSES:
        mean, sd = _mean_sd(class_deltas[name])
        print(
            "S(T0)-S(T37) %s n=%d mean=%.4f sd=%.4f"
            % (name, class_deltas[name].size, mean, sd)
        )


if __name__ == "__main__":
    main()
