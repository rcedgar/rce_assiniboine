#!/usr/bin/env python3
"""Scatter measured RNAfold Z against predicted Z for the ROD subset.

Measured Z is the Z column of RODs_subset_RNAfold_z.tsv. Predicted Z is the
Z column of RODs_predict_z_subset.tsv. Rows are joined on the label, ignoring
a leading ROD| prefix. scatter.png shows every pair. scatter_Z_-20_to_5.png
keeps pairs inside that window on both axes.
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
_DIR = _ROOT / "RODs_subset_predict_z"


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _label_key(label: str) -> str:
    return label.split("|")[-1]


def _load_z(path: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            try:
                label_col = header.index("label")
                z_col = header.index("Z")
            except ValueError:
                _die(f"ERROR: {path} is missing label or Z")
            width = max(label_col, z_col)
            for line_number, line in enumerate(handle, start=2):
                text = line.rstrip("\r\n")
                if not text:
                    continue
                fields = text.split("\t")
                if len(fields) <= width or fields[z_col] == "":
                    continue
                try:
                    value = float(fields[z_col])
                except ValueError:
                    _die(f"ERROR: {path}:{line_number} Z is not a number")
                if math.isfinite(value):
                    scores[_label_key(fields[label_col])] = value
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    if not scores:
        _die(f"ERROR: no finite Z values in {path}")
    return scores


def _write_scatter(path: Path, x: np.ndarray, y: np.ndarray, limits: tuple[float, float] | None) -> float:
    """Write the scatterplot. Return Pearson r of the plotted points."""
    if x.size < 2:
        _die(f"ERROR: need at least 2 points to plot {path.name}")
    r = float(np.corrcoef(x, y)[0, 1])
    if limits is None:
        low = float(min(x.min(), y.min()))
        high = float(max(x.max(), y.max()))
        pad = 0.05 * (high - low) if high > low else 1.0
        low -= pad
        high += pad
    else:
        low, high = limits
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=12, alpha=0.75, color="cyan", zorder=2)
    ax.axline((0, 0), slope=1, linestyle="--", color="0.35", linewidth=1, zorder=1)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Measured Z")
    ax.set_ylabel("Predicted Z")
    ax.text(
        0.04,
        0.96,
        "r = %.3f\nn = %d" % (r, x.size),
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    fig.tight_layout()
    temporary = path.with_suffix(".tmp.png")
    fig.savefig(temporary, dpi=150)
    plt.close(fig)
    temporary.replace(path)
    return r


def main() -> None:
    measured_path = _DIR / "RODs_subset_RNAfold_z.tsv"
    predicted_path = _DIR / "RODs_predict_z_subset.tsv"
    figure_path = _DIR / "scatter.png"
    measured = _load_z(measured_path)
    predicted = _load_z(predicted_path)
    labels = sorted(set(measured) & set(predicted))
    if not labels:
        _die("ERROR: no labels in both Z tables")
    only_measured = len(set(measured) - set(predicted))
    only_predicted = len(set(predicted) - set(measured))
    if only_measured or only_predicted:
        print(
            f"unmatched measured {only_measured} predicted {only_predicted}",
            file=sys.stderr,
        )
    x = np.asarray([measured[label] for label in labels], dtype=np.float64)
    y = np.asarray([predicted[label] for label in labels], dtype=np.float64)
    r = _write_scatter(figure_path, x, y, None)
    print("Pearson r = %.6g (n = %d)" % (r, x.size))
    print("wrote %s" % figure_path)
    window = (-20.0, 5.0)
    keep = (x >= window[0]) & (x <= window[1]) & (y >= window[0]) & (y <= window[1])
    window_path = _DIR / "scatter_Z_-20_to_5.png"
    r_window = _write_scatter(window_path, x[keep], y[keep], window)
    print("Z -20 to 5 Pearson r = %.6g (n = %d)" % (r_window, int(keep.sum())))
    print("wrote %s" % window_path)


if __name__ == "__main__":
    main()
