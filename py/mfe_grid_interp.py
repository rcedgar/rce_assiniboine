#!/usr/bin/env python3
"""Estimate MFE cmean and cstddev by interpolating the RNAfold grid table.

Composition uses tetrahedral interpolation on the 1/k frequency lattice.
Length is a linear blend between the tabulated lengths that bracket the query.
"""

from __future__ import annotations

import argparse
import sys
from typing import NoReturn

import numpy as np

# Offsets in (A, C, G) that partition {f in [0,1]^3 : fA + fC + fG <= 2}.
_TRUNCATED_TETS = (
    ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0)),
    ((1, 0, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1)),
    ((0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1)),
    ((0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)),
)

_ZERO = 1e-8

RowKey = tuple[int, int, int, int, int]
RowValue = tuple[int, float, float]
Support = tuple[int, int, int, int, int, int, float, float, float]


def _die(message: str) -> NoReturn:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def infer_k(freqs: np.ndarray) -> int:
    """Smallest k such that every frequency is a multiple of 1/k."""
    values = np.unique(freqs)
    for k in range(1, 65):
        scaled = k * values
        if np.all(np.abs(scaled - np.rint(scaled)) < 1e-6):
            return k
    _die("cannot infer composition grid step from frequencies")


def load_table(path: str) -> tuple[dict[RowKey, RowValue], np.ndarray, int]:
    """Load rows keyed by (L, a, c, g, t) counts. Also return sorted lengths and k."""
    try:
        data = np.loadtxt(path, comments="#", dtype=np.float64, ndmin=2)
    except OSError as exc:
        _die(f"cannot read {path}: {exc}")
    except ValueError as exc:
        _die(f"cannot parse {path}: {exc}")
    if data.size == 0 or data.shape[1] != 10:
        _die(f"expected 10 columns in {path}")

    k = infer_k(data[:, 1:5])
    rows: dict[RowKey, RowValue] = {}
    for row in data:
        length = int(round(float(row[0])))
        counts = np.rint(row[1:5] * k).astype(int)
        if np.any(np.abs(row[1:5] * k - counts) > 1e-4):
            _die(f"frequency is not a multiple of 1/{k}: {row[1:5].tolist()}")
        if int(counts.sum()) != k:
            _die(f"frequencies do not sum to 1: {row[1:5].tolist()}")
        key = (
            length,
            int(counts[0]),
            int(counts[1]),
            int(counts[2]),
            int(counts[3]),
        )
        if key in rows:
            _die(f"duplicate row L={length} counts={list(key[1:])}")
        rows[key] = (int(round(float(row[5]))), float(row[8]), float(row[9]))
    lengths = np.array(sorted({key[0] for key in rows}), dtype=np.float64)
    return rows, lengths, k


def _barycentric(point: np.ndarray, vertices: np.ndarray) -> np.ndarray | None:
    """Weights of a 3D point in a tetrahedron. vertices is (4, 3)."""
    matrix = np.empty((4, 4), dtype=np.float64)
    matrix[:3, :] = vertices.T
    matrix[3, :] = 1.0
    rhs = np.empty(4, dtype=np.float64)
    rhs[:3] = point
    rhs[3] = 1.0
    try:
        weights = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError:
        return None
    if np.max(np.abs(matrix @ weights - rhs)) > 1e-6:
        return None
    return weights


def _kuhn(frac: np.ndarray) -> list[tuple[tuple[int, int, int], float]]:
    """Kuhn tetrahedron of the unit cube containing frac."""
    order = sorted(range(3), key=lambda axis: (-float(frac[axis]), axis))
    offset = [0, 0, 0]
    offsets = [tuple(offset)]
    for axis in order:
        offset = offset.copy()
        offset[axis] += 1
        offsets.append((offset[0], offset[1], offset[2]))
    ordered = [float(frac[axis]) for axis in order]
    weights = (
        1.0 - ordered[0],
        ordered[0] - ordered[1],
        ordered[1] - ordered[2],
        ordered[2],
    )
    return list(zip(offsets, weights))


def _truncated(frac: np.ndarray, n_tets: int) -> list[tuple[tuple[int, int, int], float]]:
    """First truncated-cube tetrahedron that contains frac."""
    point = np.asarray(frac, dtype=np.float64)
    for tet in _TRUNCATED_TETS[:n_tets]:
        weights = _barycentric(point, np.asarray(tet, dtype=np.float64))
        if weights is None:
            continue
        if np.all(weights >= -_ZERO):
            return [(tet[i], float(weights[i])) for i in range(4)]
    _die("query is outside the composition simplex")


def _positive_weights(
    pairs: list[tuple[tuple[int, int, int], float]],
) -> list[tuple[tuple[int, int, int], float]]:
    kept = []
    for offset, weight in pairs:
        if weight < -1e-6:
            _die("internal error: negative composition weight")
        if weight > _ZERO:
            kept.append((offset, weight))
    total = sum(weight for _, weight in kept)
    if not kept or abs(total - 1.0) > 1e-5:
        _die("internal error: composition weights do not sum to 1")
    return [(offset, weight / total) for offset, weight in kept]


def composition_vertices(
    freqs: np.ndarray, k: int
) -> list[tuple[tuple[int, int, int, int], float]]:
    """Lattice vertices and barycentric weights for a composition on the simplex."""
    scaled = np.asarray(freqs[:3], dtype=np.float64) * k
    floors = np.floor(scaled + 1e-8).astype(int)
    frac = np.clip(scaled - floors, 0.0, None)
    tbase = int(k - int(floors.sum()))
    if tbase >= 3:
        pairs = _kuhn(frac)
    elif tbase == 2:
        pairs = _truncated(frac, 5)
    elif tbase == 1:
        pairs = _truncated(frac, 1)
    elif tbase == 0 and np.all(frac <= 1e-8):
        pairs = [((0, 0, 0), 1.0)]
    else:
        _die("query is outside the composition simplex")

    vertices = []
    for offset, weight in _positive_weights(pairs):
        counts = (
            int(floors[0] + offset[0]),
            int(floors[1] + offset[1]),
            int(floors[2] + offset[2]),
        )
        t_count = k - (counts[0] + counts[1] + counts[2])
        if min(counts[0], counts[1], counts[2], t_count) < 0:
            _die("query is outside the composition simplex")
        vertices.append(((counts[0], counts[1], counts[2], t_count), weight))
    return vertices


def _bracket(lengths: np.ndarray, length: float) -> tuple[int, int, float]:
    if length < float(lengths[0]) - 1e-6 or length > float(lengths[-1]) + 1e-6:
        _die(
            f"L={length:g} is outside tabulated range "
            f"{lengths[0]:.0f}..{lengths[-1]:.0f}"
        )
    index = int(np.searchsorted(lengths, length, side="left"))
    if index == len(lengths):
        index = len(lengths) - 1
    if abs(float(lengths[index]) - length) <= 1e-6:
        tabulated = int(lengths[index])
        return tabulated, tabulated, 0.0
    if index == 0:
        _die(
            f"L={length:g} is outside tabulated range "
            f"{lengths[0]:.0f}..{lengths[-1]:.0f}"
        )
    lo = int(lengths[index - 1])
    hi = int(lengths[index])
    return lo, hi, (length - lo) / (hi - lo)


def _missing_message(k: int, keys: list[RowKey]) -> str:
    parts = []
    for length, a_count, c_count, g_count, t_count in keys:
        parts.append(
            "L=%d fA=%.3f fC=%.3f fG=%.3f fT=%.3f"
            % (
                length,
                a_count / k,
                c_count / k,
                g_count / k,
                t_count / k,
            )
        )
    return "missing grid row " + "; ".join(parts)


def interpolate(
    rows: dict[RowKey, RowValue],
    lengths: np.ndarray,
    k: int,
    length: float,
    freqs: np.ndarray,
) -> tuple[float, float, list[Support]]:
    """Return cmean, cstddev, and the weighted source rows."""
    vertices = composition_vertices(freqs, k)
    lo, hi, t_length = _bracket(lengths, length)
    if lo == hi:
        slices = [(lo, 1.0)]
    else:
        slices = [(lo, 1.0 - t_length), (hi, t_length)]

    missing: list[RowKey] = []
    support: list[Support] = []
    for tabulated, length_weight in slices:
        for counts, composition_weight in vertices:
            weight = length_weight * composition_weight
            if weight <= _ZERO:
                continue
            key = (tabulated, counts[0], counts[1], counts[2], counts[3])
            row = rows.get(key)
            if row is None:
                missing.append(key)
                continue
            n_samples, cmean, cstddev = row
            support.append(
                (
                    tabulated,
                    counts[0],
                    counts[1],
                    counts[2],
                    counts[3],
                    n_samples,
                    cmean,
                    cstddev,
                    weight,
                )
            )
    if missing:
        _die(_missing_message(k, missing))

    support.sort(key=lambda item: item[:5])
    weights = np.array([item[8] for item in support], dtype=np.float64)
    total = float(weights.sum())
    if abs(total - 1.0) > 1e-5:
        _die("internal error: weights do not sum to 1")
    weights /= total
    cmeans = np.array([item[6] for item in support], dtype=np.float64)
    cstddevs = np.array([item[7] for item in support], dtype=np.float64)
    normalized = [
        (
            item[0],
            item[1],
            item[2],
            item[3],
            item[4],
            item[5],
            item[6],
            item[7],
            float(weights[index]),
        )
        for index, item in enumerate(support)
    ]
    return float(np.dot(weights, cmeans)), float(np.dot(weights, cstddevs)), normalized


def _normalize_freqs(fA: float, fC: float, fG: float, fT: float) -> np.ndarray:
    freqs = np.array([fA, fC, fG, fT], dtype=np.float64)
    if not np.all(np.isfinite(freqs)) or np.any(freqs < 0.0):
        _die("frequencies must be finite and >= 0")
    total = float(freqs.sum())
    if total < 0.99 or total > 1.01:
        _die(f"frequencies sum to {total}, expected 0.99..1.01")
    return freqs / total


def _format_length(length: float) -> str:
    rounded = round(length)
    if abs(length - rounded) <= 1e-6:
        return str(int(rounded))
    return f"{length:.6g}"


def _format_freq(value: float) -> str:
    return f"{value:.6g}"


def _print_report(
    length: float,
    freqs: np.ndarray,
    k: int,
    cmean: float,
    cstddev: float,
    support: list[Support],
) -> None:
    sys.stdout.write("L\tfA\tfC\tfG\tfT\tn\tcmean\tcstddev\tweight\n")
    for item in support:
        tabulated, a_count, c_count, g_count, t_count, n_samples, row_cmean, row_csd, weight = item
        sys.stdout.write(
            "%d\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3g\t%.3g\t%.4f\n"
            % (
                tabulated,
                a_count / k,
                c_count / k,
                g_count / k,
                t_count / k,
                n_samples,
                row_cmean,
                row_csd,
                weight,
            )
        )
    sys.stdout.write(
        "%s\t%s\t%s\t%s\t%s\t\t%.3g\t%.3g\t\n"
        % (
            _format_length(length),
            _format_freq(float(freqs[0])),
            _format_freq(float(freqs[1])),
            _format_freq(float(freqs[2])),
            _format_freq(float(freqs[3])),
            cmean,
            cstddev,
        )
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Estimate MFE cmean and cstddev by interpolating the RNAfold grid table."
    )
    parser.add_argument("L", type=float, help="Sequence length")
    parser.add_argument("fA", type=float, help="Frequency of A")
    parser.add_argument("fC", type=float, help="Frequency of C")
    parser.add_argument("fG", type=float, help="Frequency of G")
    parser.add_argument("fT", type=float, help="Frequency of T")
    parser.add_argument(
        "--table",
        default="../output/rnafold_grid_mean_var.tsv",
        help="Grid table (default: ../output/rnafold_grid_mean_var.tsv)",
    )
    args = parser.parse_args(argv)

    if not np.isfinite(args.L) or args.L <= 0.0:
        _die("L must be finite and > 0")
    freqs = _normalize_freqs(args.fA, args.fC, args.fG, args.fT)
    rows, lengths, k = load_table(args.table)
    cmean, cstddev, support = interpolate(rows, lengths, k, args.L, freqs)
    _print_report(args.L, freqs, k, cmean, cstddev, support)


if __name__ == "__main__":
    main()
