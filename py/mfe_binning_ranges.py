#!/usr/bin/env python3
"""Half-open L and base-count ranges for each leaf of an MFE binning tree.

The tree splits on log(L), GC fraction, AT skew, and GC skew. Each leaf is the
product of half-open intervals inherited from those splits. Count ranges are
the infimum and supremum of NA, NC, NG, and NT over that cell, with counts
real-valued, nonnegative, and summing to L.
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path
from typing import NamedTuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import mfe_binning_features as feat

OUTPUT_HEADER = (
    "bin\tL_lo\tL_hi\tNA_lo\tNA_hi\tNC_lo\tNC_hi\tNG_lo\tNG_hi\tNT_lo\tNT_hi"
)
_VERTEX_TOL = 1e-8


class ParameterRanges(NamedTuple):
    """Infimum and supremum of each parameter on one leaf. The interval is [lo, hi)."""

    L_lo: float
    L_hi: float
    NA_lo: float
    NA_hi: float
    NC_lo: float
    NC_hi: float
    NG_lo: float
    NG_hi: float
    NT_lo: float
    NT_hi: float


def all_feature_cells(tree: feat.Tree) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(lo, hi)`` with shape ``(n_nodes, 4)``.

    Feature ``k`` on a node lies in ``[lo, hi)``. The root starts unbounded.
    A left child tightens the exclusive upper bound to the parent's threshold.
    A right child tightens the inclusive lower bound to that threshold.
    """
    n_nodes = int(tree.is_leaf.shape[0])
    lo = np.full((n_nodes, feat.N_FEATURES), -np.inf, dtype=np.float64)
    hi = np.full((n_nodes, feat.N_FEATURES), np.inf, dtype=np.float64)
    stack = [0]
    while stack:
        node = stack.pop()
        if tree.is_leaf[node]:
            continue
        axis = int(tree.axis[node])
        thr = float(tree.thr[node])
        left = int(tree.left[node])
        right = int(tree.right[node])
        lo[left] = lo[node]
        hi[left] = hi[node]
        hi[left, axis] = min(hi[left, axis], thr)
        lo[right] = lo[node]
        hi[right] = hi[node]
        lo[right, axis] = max(lo[right, axis], thr)
        stack.append(right)
        stack.append(left)
    return lo, hi


def _composition_vertices(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Vertices of the closed composition box cut by the count simplex.

    ``lo`` and ``hi`` are the x1, x2, x3 bounds. The simplex requires
    ``x1 in [0, 1]``, ``|x2| <= 1 - x1``, and ``|x3| <= x1``.
    """
    box_lo = np.array(
        [max(lo[0], 0.0), max(lo[1], -1.0), max(lo[2], -1.0)],
        dtype=np.float64,
    )
    box_hi = np.array(
        [min(hi[0], 1.0), min(hi[1], 1.0), min(hi[2], 1.0)],
        dtype=np.float64,
    )
    if np.any(box_lo > box_hi + _VERTEX_TOL):
        raise ValueError("feature cell does not meet the composition domain")

    # Each entry is (a, b) for the inequality a · (x1, x2, x3) <= b.
    inequalities: list[tuple[np.ndarray, float]] = [
        (np.array([1.0, 0.0, 0.0]), float(box_hi[0])),
        (np.array([-1.0, 0.0, 0.0]), float(-box_lo[0])),
        (np.array([0.0, 1.0, 0.0]), float(box_hi[1])),
        (np.array([0.0, -1.0, 0.0]), float(-box_lo[1])),
        (np.array([0.0, 0.0, 1.0]), float(box_hi[2])),
        (np.array([0.0, 0.0, -1.0]), float(-box_lo[2])),
        (np.array([1.0, 1.0, 0.0]), 1.0),   # x2 <= 1 - x1
        (np.array([1.0, -1.0, 0.0]), 1.0),  # x2 >= x1 - 1
        (np.array([-1.0, 0.0, 1.0]), 0.0),  # x3 <= x1
        (np.array([-1.0, 0.0, -1.0]), 0.0),  # x3 >= -x1
    ]
    vertices: list[np.ndarray] = []
    for chosen in combinations(range(len(inequalities)), 3):
        matrix = np.stack([inequalities[i][0] for i in chosen])
        rhs = np.array([inequalities[i][1] for i in chosen], dtype=np.float64)
        det = float(np.linalg.det(matrix))
        if not np.isfinite(det) or abs(det) < 1e-12:
            continue
        point = np.linalg.solve(matrix, rhs)
        if not np.all(np.isfinite(point)):
            continue
        if any(float(a @ point) > bound + _VERTEX_TOL for a, bound in inequalities):
            continue
        if any(float(np.linalg.norm(point - kept)) < 1e-8 for kept in vertices):
            continue
        vertices.append(point)
    if not vertices:
        raise ValueError("composition cell has no vertices")
    return np.stack(vertices)


def _snap_integer(value: float) -> float:
    """Replace a value that is only a rounding error away from an integer."""
    if not np.isfinite(value):
        return float(value)
    nearest = float(np.rint(value))
    if abs(value - nearest) <= 1e-8 * max(1.0, abs(value)):
        return nearest
    return float(value)


def integer_half_open(lo: float, hi: float) -> tuple[float, float]:
    """Integers ``n`` with ``lo <= n < hi``, as ``[ceil(lo), ceil(hi))``.

    ``inf`` stays ``inf``. A value within a few ulps of an integer is snapped
    first so ``ceil(20.000000000000004)`` does not become 21.
    """
    if np.isfinite(lo):
        lo_i = float(np.ceil(_snap_integer(lo)))
    else:
        lo_i = float("-inf")
    if not np.isfinite(hi):
        return lo_i, float("inf") if hi > 0 else float("-inf")
    return lo_i, float(np.ceil(_snap_integer(hi)))


def _count_span(fraction: np.ndarray, length_lo: float, length_hi: float) -> tuple[float, float]:
    """Inf and sup of ``fraction * L`` for L in ``[length_lo, length_hi)`` and ``fraction >= 0``."""
    low = float(np.min(fraction))
    high = float(np.max(fraction))
    if low < 0.0 and low > -1e-9:
        low = 0.0
    if high < 0.0 and high > -1e-9:
        high = 0.0
    count_lo = low * length_lo
    if np.isfinite(length_hi):
        count_hi = high * length_hi
    elif high > 0.0:
        count_hi = float("inf")
    else:
        count_hi = 0.0
    return count_lo, count_hi


def parameter_ranges(lo: np.ndarray, hi: np.ndarray) -> ParameterRanges:
    """Map one feature cell ``[lo, hi)`` to integer L and base-count ranges.

    ``lo`` and ``hi`` have length 4: log(L), GC fraction, AT skew, GC skew.
    Returned bounds are integers, or ``inf``. ``[lo, hi)`` contains every
    integer value of that parameter inside the real cell.
    """
    lo = np.asarray(lo, dtype=np.float64).reshape(4)
    hi = np.asarray(hi, dtype=np.float64).reshape(4)
    length_lo = float(np.exp(lo[0])) if np.isfinite(lo[0]) else 0.0
    length_hi = float(np.exp(hi[0])) if np.isfinite(hi[0]) else float("inf")
    vertices = _composition_vertices(lo[1:], hi[1:])
    x1 = vertices[:, 0]
    x2 = vertices[:, 1]
    x3 = vertices[:, 2]
    fractions = (
        (1.0 - x1 + x2) / 2.0,  # NA / L
        (x1 - x3) / 2.0,        # NC / L
        (x1 + x3) / 2.0,        # NG / L
        (1.0 - x1 - x2) / 2.0,  # NT / L
    )
    spans = [_count_span(fraction, length_lo, length_hi) for fraction in fractions]
    ends = [integer_half_open(length_lo, length_hi)]
    ends.extend(integer_half_open(span_lo, span_hi) for span_lo, span_hi in spans)
    return ParameterRanges(*(value for pair in ends for value in pair))


def leaf_ranges(tree: feat.Tree) -> list[ParameterRanges]:
    """One range row per leaf, in increasing node_id. Row ``i`` is bin ``i``."""
    lo, hi = all_feature_cells(tree)
    leaves = np.flatnonzero(tree.is_leaf)
    return [parameter_ranges(lo[int(leaf)], hi[int(leaf)]) for leaf in leaves]


def _fmt(value: float) -> str:
    if not np.isfinite(value):
        return "inf" if value > 0 else "-inf"
    return str(int(value))


def write_ranges(path: str, rows: list[ParameterRanges]) -> None:
    """Write the range TSV. ``path`` of ``-`` writes to stdout."""
    own_handle = path != "-"
    handle = open(path, "w", encoding="utf-8", newline="\n") if own_handle else sys.stdout
    try:
        handle.write(OUTPUT_HEADER + "\n")
        for bin_index, row in enumerate(rows):
            fields = [str(bin_index)]
            fields.extend(_fmt(value) for value in row)
            handle.write("\t".join(fields) + "\n")
    finally:
        if own_handle:
            handle.close()


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Write [lo, hi) ranges for each MFE binning leaf.")
    parser.add_argument("--tree", required=True, help="Tree TSV from mfe_binning_build_bins.py")
    parser.add_argument("--output", default="-", help="Output TSV, or - for stdout")
    args = parser.parse_args(argv)
    try:
        tree = feat.load_tree(args.tree)
        rows = leaf_ranges(tree)
    except ValueError as exc:
        _die(str(exc))
    write_ranges(args.output, rows)


if __name__ == "__main__":
    main()
