#!/usr/bin/env python3
"""Build an adaptive (length, composition) partition tree for MFE records.

The tree is grown on a training split only, with an explicit stack (no
recursion). Leaf statistics are then filled for the holdout by routing those
rows through :func:`mfe_binning_features.assign`.

Variance is the sample variance (ddof=1). Skew and excess kurtosis use
bias-uncorrected central moments (divide by n, not n-1). The split's residual
variance is SSE / (n - 5), one degree of freedom for the intercept and one
for each of the four slopes.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import NamedTuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import mfe_binning_features as feat

CHUNK = 500_000
HEADER_TOKENS = {"MFE", "L", "length", "NA", "NC", "NG", "NT"}
_COUNT_ATOL = 1e-4


class SplitDecision(NamedTuple):
    """Result of :func:`evaluate_node`.

    ``stop_reason`` is set only for leaves. ``cut`` is the left-child size
    when a split is chosen. ``bias_var`` and ``s2`` are NaN when the node
    stopped before the regression (min size or max depth).
    """

    is_leaf: bool
    stop_reason: str | None
    axis: int
    threshold: float
    cut: int
    bias_var: float
    s2: float
    slopes: np.ndarray | None


def _chunks(idx: np.ndarray, chunk: int = CHUNK):
    n = int(idx.shape[0])
    for start in range(0, n, chunk):
        yield idx[start : start + chunk]


def stable_threshold(lo, hi) -> np.float64:
    """Return a float64 threshold with ``lo < thr <= hi`` that survives ``%.17g``.

    The value returned is the number obtained by formatting and parsing, so the
    in-memory cut and the reloaded TSV cut are the same bits.
    """
    lo = np.float64(lo)
    hi = np.float64(hi)
    if not (lo < hi):
        raise ValueError("threshold bounds must satisfy lo < hi")
    thr = (lo + hi) * np.float64(0.5)
    if not (lo < thr <= hi):
        thr = np.nextafter(lo, np.float64(np.inf))
    for _ in range(64):
        if not (lo < thr <= hi):
            break
        parsed = np.float64(format(float(thr), ".17g"))
        if lo < parsed <= hi:
            return parsed
        if parsed <= lo:
            thr = np.nextafter(lo, np.float64(np.inf))
        else:
            thr = np.nextafter(thr, lo)
    parsed = np.float64(format(float(hi), ".17g"))
    if lo < parsed <= hi:
        return parsed
    thr = np.nextafter(lo, np.float64(np.inf))
    steps = 0
    while thr <= hi and steps < 10000:
        parsed = np.float64(format(float(thr), ".17g"))
        if lo < parsed <= hi:
            return parsed
        nxt = np.nextafter(thr, np.float64(np.inf))
        if nxt == thr:
            break
        thr = nxt
        steps += 1
    raise ValueError("no threshold in (lo, hi] round-trips through %.17g")


def find_cut(values: np.ndarray, n_min: int):
    """Nearest-to-median cut with both sides at least ``n_min`` and a strict gap.

    Returns ``(cut, threshold)`` where ``cut`` is the left-child size on the
    sorted axis, or None when every gap is tied or too close to an edge.
    """
    values = np.asarray(values, dtype=np.float64)
    n = int(values.shape[0])
    if n < 2 * n_min:
        return None
    order = np.argsort(values, kind="mergesort")
    sorted_vals = values[order]
    cuts = np.flatnonzero(sorted_vals[1:] > sorted_vals[:-1]) + 1
    cuts = cuts[(cuts >= n_min) & (cuts <= n - n_min)]
    if cuts.size == 0:
        return None
    target = n // 2
    best = int(np.argmin(np.abs(cuts.astype(np.int64) - target)))
    cut = int(cuts[best])
    lo = np.float64(sorted_vals[cut - 1])
    hi = np.float64(sorted_vals[cut])
    return cut, stable_threshold(lo, hi)


def ols_bias(X: np.ndarray, y: np.ndarray, idx: np.ndarray):
    """Centered OLS of ``y`` on the four features, via the 4x4 normal matrix.

    Returns ``(slopes, explained, bias_var, s2)``. ``explained[i]`` is
    ``b_i^2 * Var(x_i)`` with ddof=1. ``s2`` is ``SSE / (n - 5)``. Rows are
    visited in chunks so ``X`` itself is not copied.
    """
    n = int(idx.shape[0])
    if n <= 5:
        raise ValueError("OLS needs more than 5 rows")
    sum_x = np.zeros(feat.N_FEATURES, dtype=np.float64)
    sum_y = np.float64(0.0)
    for sl in _chunks(idx):
        block = X[sl]
        sum_x += block.sum(axis=0, dtype=np.float64)
        sum_y += y[sl].astype(np.float64, copy=False).sum(dtype=np.float64)
    mean_x = sum_x / n
    mean_y = sum_y / n
    gram = np.zeros((feat.N_FEATURES, feat.N_FEATURES), dtype=np.float64)
    xty = np.zeros(feat.N_FEATURES, dtype=np.float64)
    for sl in _chunks(idx):
        centered_x = X[sl] - mean_x
        centered_y = y[sl].astype(np.float64, copy=False) - mean_y
        gram += centered_x.T @ centered_x
        xty += centered_x.T @ centered_y
    slopes, _, _, _ = np.linalg.lstsq(gram, xty, rcond=None)
    slopes = np.asarray(slopes, dtype=np.float64).reshape(feat.N_FEATURES)
    sse = np.float64(0.0)
    for sl in _chunks(idx):
        centered_x = X[sl] - mean_x
        centered_y = y[sl].astype(np.float64, copy=False) - mean_y
        resid = centered_y - centered_x @ slopes
        sse += resid @ resid
    if sse < 0.0:
        sse = np.float64(0.0)
    var_x = np.diag(gram) / (n - 1)
    explained = slopes * slopes * var_x
    bias_var = float(np.sum(explained))
    s2 = float(sse / (n - 5))
    return slopes, explained, bias_var, s2


def _leaf_decision(reason: str, bias_var: float = np.nan, s2: float = np.nan, slopes=None):
    return SplitDecision(
        is_leaf=True,
        stop_reason=reason,
        axis=-1,
        threshold=float("nan"),
        cut=0,
        bias_var=float(bias_var),
        s2=float(s2),
        slopes=slopes,
    )


def evaluate_node(X, y, idx, depth: int, n_min: int, tau: float, max_depth: int) -> SplitDecision:
    """Decide whether the rows ``idx`` are a leaf, or return the split.

    Stop order is ``min_n``, then ``max_depth``, then a regression. A leaf is
    ``bias_ok`` when ``bias_var / s2 < tau``. Otherwise axes are tried from
    the largest ``b_i^2 Var(x_i)`` downward. ``no_split`` means no axis had a
    cut that keeps both children at least ``n_min``.
    """
    idx = np.asarray(idx)
    n = int(idx.shape[0])
    if n < 2 * n_min:
        return _leaf_decision("min_n")
    if depth >= max_depth:
        return _leaf_decision("max_depth")
    if n <= 5:
        return _leaf_decision("no_split")
    slopes, explained, bias_var, s2 = ols_bias(X, y, idx)
    if s2 > 0.0 and (bias_var / s2) < tau:
        return _leaf_decision("bias_ok", bias_var, s2, slopes)
    axes = np.argsort(-explained, kind="mergesort")
    for axis in axes:
        axis = int(axis)
        found = find_cut(X[idx, axis], n_min)
        if found is None:
            continue
        cut, thr = found
        return SplitDecision(
            is_leaf=False,
            stop_reason=None,
            axis=axis,
            threshold=float(thr),
            cut=int(cut),
            bias_var=bias_var,
            s2=s2,
            slopes=slopes,
        )
    return _leaf_decision("no_split", bias_var, s2, slopes)


def moment_stats(y: np.ndarray, *, with_quantiles: bool = True) -> dict:
    """Mean, sample variance, and bias-uncorrected skew and excess kurtosis.

    Skew and excess kurtosis are blank (NaN) when ``n < 2`` or the second
    central moment is 0. Quantiles use numpy's linear method.
    """
    values = np.asarray(y, dtype=np.float64).reshape(-1)
    n = int(values.shape[0])
    nan = float("nan")
    quantiles = np.full(4, nan, dtype=np.float64)
    if n == 0:
        return {
            "n": 0,
            "mean": nan,
            "var": nan,
            "skew": nan,
            "exkurt": nan,
            "quantiles": quantiles,
        }
    mean = float(values.mean())
    if with_quantiles:
        quantiles = np.quantile(values, [0.001, 0.01, 0.99, 0.999]).astype(np.float64)
    if n < 2:
        return {
            "n": n,
            "mean": mean,
            "var": nan,
            "skew": nan,
            "exkurt": nan,
            "quantiles": quantiles,
        }
    centered = values - mean
    m2 = float(np.dot(centered, centered) / n)
    var = float(np.dot(centered, centered) / (n - 1))
    if m2 == 0.0:
        skew = nan
        exkurt = nan
    else:
        m3 = float(np.dot(centered * centered, centered) / n)
        m4 = float(np.dot(centered * centered, centered * centered) / n)
        skew = m3 / (m2 ** 1.5)
        exkurt = m4 / (m2 * m2) - 3.0
    return {
        "n": n,
        "mean": mean,
        "var": var,
        "skew": float(skew),
        "exkurt": float(exkurt),
        "quantiles": quantiles,
    }


def _bbox(X: np.ndarray, idx: np.ndarray):
    lo = np.empty(feat.N_FEATURES, dtype=np.float64)
    hi = np.empty(feat.N_FEATURES, dtype=np.float64)
    for axis in range(feat.N_FEATURES):
        column = X[idx, axis]
        lo[axis] = column.min()
        hi[axis] = column.max()
    return lo, hi


def build_tree(X, y, n_min: int, tau: float, max_depth: int) -> feat.Tree:
    """Partition training rows and return a tree with training leaf statistics.

    Holdout columns are left empty (NaN, and ``n_hold`` 0 on leaves) for
    :func:`fill_holdout_stats` to fill. ``X`` is not copied; only the index
    array is reordered.
    """
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y)
    if X.ndim != 2 or X.shape[1] != feat.N_FEATURES:
        raise ValueError("X must have shape (n, 4)")
    n = int(X.shape[0])
    if y.shape[0] != n:
        raise ValueError("X and y have different lengths")
    if n == 0:
        raise ValueError("training set is empty")
    if n_min < 1:
        raise ValueError("n_min must be >= 1")
    if max_depth < 0:
        raise ValueError("max_depth must be >= 0")

    idx_dtype = np.int32 if n <= np.iinfo(np.int32).max else np.int64
    order = np.arange(n, dtype=idx_dtype)
    # (node_id, start, end, depth, parent). Push right then left so left is next.
    stack: list[tuple[int, int, int, int, int]] = [(0, 0, n, 0, -1)]
    next_id = 1
    nodes: dict[int, dict] = {}
    processed = 0
    n_leaves = 0

    while stack:
        node_id, start, end, depth, parent = stack.pop()
        idx = order[start:end]
        decision = evaluate_node(X, y, idx, depth, n_min, tau, max_depth)
        processed += 1
        record = {
            "parent": parent,
            "depth": depth,
            "is_leaf": decision.is_leaf,
            "axis": -1,
            "thr": float("nan"),
            "left": -1,
            "right": -1,
            "stop_reason": "",
            "stats": None,
            "bbox_min": None,
            "bbox_max": None,
        }
        if decision.is_leaf:
            n_leaves += 1
            record["stop_reason"] = decision.stop_reason or ""
            record["stats"] = moment_stats(y[idx])
            record["bbox_min"], record["bbox_max"] = _bbox(X, idx)
        else:
            taken = idx.copy()
            values = X[taken, decision.axis]
            left_mask = values < decision.threshold
            n_left = int(np.count_nonzero(left_mask))
            n_right = int(taken.shape[0] - n_left)
            if n_left != decision.cut or n_left < n_min or n_right < n_min:
                raise RuntimeError(
                    f"split at node {node_id} produced children {n_left} and {n_right}"
                )
            order[start : start + n_left] = taken[left_mask]
            order[start + n_left : end] = taken[~left_mask]
            left_id = next_id
            right_id = next_id + 1
            next_id += 2
            record["axis"] = decision.axis
            record["thr"] = decision.threshold
            record["left"] = left_id
            record["right"] = right_id
            stack.append((right_id, start + n_left, end, depth + 1, node_id))
            stack.append((left_id, start, start + n_left, depth + 1, node_id))
        nodes[node_id] = record
        if processed % 100 == 0:
            print(
                f"nodes={processed} leaves={n_leaves} queue={len(stack)}",
                file=sys.stderr,
                flush=True,
            )

    n_nodes = next_id
    if len(nodes) != n_nodes or sorted(nodes) != list(range(n_nodes)):
        raise RuntimeError("node ids are not the contiguous range 0 .. n_nodes-1")
    return _tree_from_records(nodes, n_nodes)


def _nan_array(n: int) -> np.ndarray:
    return np.full(n, np.nan, dtype=np.float64)


def _tree_from_records(nodes: dict[int, dict], n_nodes: int) -> feat.Tree:
    axis = np.full(n_nodes, -1, dtype=np.int32)
    thr = _nan_array(n_nodes)
    left = np.full(n_nodes, -1, dtype=np.int32)
    right = np.full(n_nodes, -1, dtype=np.int32)
    parent = np.full(n_nodes, -1, dtype=np.int32)
    is_leaf = np.zeros(n_nodes, dtype=bool)
    depth = np.zeros(n_nodes, dtype=np.int32)
    stop_reason = np.empty(n_nodes, dtype=object)
    n_train = _nan_array(n_nodes)
    mean_train = _nan_array(n_nodes)
    var_train = _nan_array(n_nodes)
    skew_train = _nan_array(n_nodes)
    exkurt_train = _nan_array(n_nodes)
    q001 = _nan_array(n_nodes)
    q01 = _nan_array(n_nodes)
    q99 = _nan_array(n_nodes)
    q999 = _nan_array(n_nodes)
    n_hold = _nan_array(n_nodes)
    bbox_min = np.full((n_nodes, feat.N_FEATURES), np.nan, dtype=np.float64)
    bbox_max = np.full((n_nodes, feat.N_FEATURES), np.nan, dtype=np.float64)
    for node_id, record in nodes.items():
        axis[node_id] = record["axis"]
        thr[node_id] = record["thr"]
        left[node_id] = record["left"]
        right[node_id] = record["right"]
        parent[node_id] = record["parent"]
        is_leaf[node_id] = record["is_leaf"]
        depth[node_id] = record["depth"]
        stop_reason[node_id] = record["stop_reason"]
        if record["is_leaf"]:
            stats = record["stats"]
            n_train[node_id] = stats["n"]
            mean_train[node_id] = stats["mean"]
            var_train[node_id] = stats["var"]
            skew_train[node_id] = stats["skew"]
            exkurt_train[node_id] = stats["exkurt"]
            q001[node_id] = stats["quantiles"][0]
            q01[node_id] = stats["quantiles"][1]
            q99[node_id] = stats["quantiles"][2]
            q999[node_id] = stats["quantiles"][3]
            n_hold[node_id] = 0.0
            bbox_min[node_id] = record["bbox_min"]
            bbox_max[node_id] = record["bbox_max"]
    return feat.Tree(
        axis=axis,
        thr=thr,
        left=left,
        right=right,
        parent=parent,
        is_leaf=is_leaf,
        depth=depth,
        stop_reason=stop_reason,
        n_train=n_train,
        mean_train=mean_train,
        var_train=var_train,
        skew_train=skew_train,
        exkurt_train=exkurt_train,
        q001_train=q001,
        q01_train=q01,
        q99_train=q99,
        q999_train=q999,
        n_hold=n_hold,
        mean_hold=_nan_array(n_nodes),
        var_hold=_nan_array(n_nodes),
        skew_hold=_nan_array(n_nodes),
        exkurt_hold=_nan_array(n_nodes),
        bbox_min=bbox_min,
        bbox_max=bbox_max,
    )


def fill_holdout_stats(tree: feat.Tree, X_hold, y_hold) -> None:
    """Route holdout rows with :func:`feat.assign` and store per-leaf moments."""
    X_hold = np.asarray(X_hold, dtype=np.float64)
    y_hold = np.asarray(y_hold)
    n = int(X_hold.shape[0])
    leaves = tree.is_leaf
    tree.n_hold[:] = np.nan
    tree.n_hold[leaves] = 0.0
    tree.mean_hold[:] = np.nan
    tree.var_hold[:] = np.nan
    tree.skew_hold[:] = np.nan
    tree.exkurt_hold[:] = np.nan
    if n == 0:
        return
    leaf_ids = feat.assign(tree, X_hold)
    order = np.argsort(leaf_ids, kind="mergesort")
    sorted_ids = leaf_ids[order]
    sorted_y = y_hold[order]
    change = np.flatnonzero(sorted_ids[1:] != sorted_ids[:-1]) + 1
    starts = np.concatenate((np.array([0], dtype=np.int64), change.astype(np.int64)))
    ends = np.concatenate((change.astype(np.int64), np.array([n], dtype=np.int64)))
    for start, end in zip(starts, ends):
        leaf_id = int(sorted_ids[start])
        stats = moment_stats(sorted_y[start:end], with_quantiles=False)
        tree.n_hold[leaf_id] = stats["n"]
        tree.mean_hold[leaf_id] = stats["mean"]
        tree.var_hold[leaf_id] = stats["var"]
        tree.skew_hold[leaf_id] = stats["skew"]
        tree.exkurt_hold[leaf_id] = stats["exkurt"]


def holdout_z_moments(tree: feat.Tree, X_hold, y_hold) -> dict | None:
    """Moments of ``(MFE - mean_train) / sd_train`` on holdout rows.

    Rows whose leaf has non-positive or non-finite training variance are
    skipped. Returns None when fewer than two rows remain. ``var`` is ddof=1;
    skew and excess kurtosis match :func:`moment_stats`.
    """
    y_hold = np.asarray(y_hold, dtype=np.float64).reshape(-1)
    if y_hold.shape[0] == 0:
        return None
    leaf_ids = feat.assign(tree, np.asarray(X_hold, dtype=np.float64))
    mean = tree.mean_train[leaf_ids]
    var = tree.var_train[leaf_ids]
    ok = np.isfinite(y_hold) & np.isfinite(mean) & np.isfinite(var) & (var > 0.0)
    if int(np.count_nonzero(ok)) < 2:
        return None
    z = (y_hold[ok] - mean[ok]) / np.sqrt(var[ok])
    return moment_stats(z, with_quantiles=False)


def format_summary(tree: feat.Tree, z_stats: dict | None) -> str:
    """Text summary of leaf sizes, stop reasons, and holdout z moments."""
    leaves = tree.is_leaf
    sizes = tree.n_train[leaves]
    depth = int(np.max(tree.depth)) if tree.depth.size else 0
    n_leaves = int(np.count_nonzero(leaves))
    reasons = tree.stop_reason[leaves]
    reason_counts = [
        (name, int(np.sum(reasons == name)))
        for name in ("bias_ok", "min_n", "max_depth", "no_split")
    ]
    lines = [
        f"leaves={n_leaves} depth={depth}",
        (
            f"leaf_n_train min={np.min(sizes):.0f} "
            f"median={np.median(sizes):.0f} max={np.max(sizes):.0f}"
        ),
        "stop_reason " + " ".join(f"{name}={count}" for name, count in reason_counts),
    ]
    if z_stats is None:
        lines.append("holdout_z mean=NA sd=NA skew=NA exkurt=NA")
    else:
        sd = float(np.sqrt(z_stats["var"])) if np.isfinite(z_stats["var"]) else float("nan")
        lines.append(
            f"holdout_z mean={z_stats['mean']:.6g} sd={sd:.6g} "
            f"skew={z_stats['skew']:.6g} exkurt={z_stats['exkurt']:.6g}"
        )
    return "\n".join(lines)


def _fmt_float(value: float) -> str:
    if not np.isfinite(value):
        return ""
    return "%.17g" % value


def _fmt_count(value: float) -> str:
    if not np.isfinite(value):
        return ""
    return str(int(round(float(value))))


def write_tree_tsv(path: str, tree: feat.Tree) -> None:
    """Write one row per node, sorted by node_id. Leaf stats are blank on internal nodes."""
    n_nodes = int(tree.is_leaf.shape[0])
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(feat.TREE_COLUMNS) + "\n")
        for node_id in range(n_nodes):
            leaf = bool(tree.is_leaf[node_id])
            fields = [
                str(node_id),
                str(int(tree.parent[node_id])),
                str(int(tree.left[node_id])),
                str(int(tree.right[node_id])),
                str(int(tree.axis[node_id])),
                "" if leaf else ("%.17g" % tree.thr[node_id]),
                str(int(tree.depth[node_id])),
                "1" if leaf else "0",
                tree.stop_reason[node_id] if leaf else "",
                _fmt_count(tree.n_train[node_id]),
                _fmt_float(tree.mean_train[node_id]),
                _fmt_float(tree.var_train[node_id]),
                _fmt_float(tree.skew_train[node_id]),
                _fmt_float(tree.exkurt_train[node_id]),
                _fmt_float(tree.q001_train[node_id]),
                _fmt_float(tree.q01_train[node_id]),
                _fmt_float(tree.q99_train[node_id]),
                _fmt_float(tree.q999_train[node_id]),
                _fmt_count(tree.n_hold[node_id]),
                _fmt_float(tree.mean_hold[node_id]),
                _fmt_float(tree.var_hold[node_id]),
                _fmt_float(tree.skew_hold[node_id]),
                _fmt_float(tree.exkurt_hold[node_id]),
            ]
            for axis in range(feat.N_FEATURES):
                fields.append(_fmt_float(tree.bbox_min[node_id, axis]))
                fields.append(_fmt_float(tree.bbox_max[node_id, axis]))
            handle.write("\t".join(fields) + "\n")


def write_meta(path: str, *, parameters: dict, rows: dict, tree: feat.Tree) -> None:
    """Write the JSON sidecar describing features, parameters, and row counts."""
    payload = {
        "features": [
            {"name": name, "definition": definition}
            for name, definition in feat.FEATURE_DEFINITIONS
        ],
        "feature_dtype": "float64",
        "parameters": parameters,
        "rows": rows,
        "tree_depth": int(np.max(tree.depth)) if tree.depth.size else 0,
        "n_leaves": int(np.count_nonzero(tree.is_leaf)),
        "n_nodes": int(tree.is_leaf.shape[0]),
        "split_convention": feat.SPLIT_CONVENTION,
        "variance": "ddof=1",
        "skew_kurtosis": "bias-uncorrected central moments (divide by n)",
        "residual_variance": "SSE / (n - 5)",
        "quantiles": "linear quantiles at 0.001, 0.01, 0.99, 0.999",
    }
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def holdout_split(n: int, holdout_frac: float, seed: int):
    """Permutation split. ``n_hold = round(holdout_frac * n)``."""
    if not 0.0 <= float(holdout_frac) <= 1.0:
        raise ValueError("holdout_frac must be between 0 and 1")
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_hold = int(round(float(holdout_frac) * n))
    return perm[n_hold:], perm[:n_hold]


def _to_counts(values: np.ndarray):
    """Round a float column to int64. The mask is True where that is not exact."""
    finite = np.isfinite(values)
    safe = np.where(finite, values, 0.0)
    rounded = np.rint(safe)
    inexact = (~finite) | (np.abs(safe - rounded) > _COUNT_ATOL)
    return rounded.astype(np.int64), inexact


def filter_rows(L, MFE, NA, NC, NG, NT, length2=None):
    """Drop non-positive L, count sums other than L, non-finite MFE, bad length2.

    Reason counts can overlap; ``dropped`` counts each bad row once. Non-integral
    count fields are included in ``sum!=L``.
    """
    length = np.asarray(L)
    mfe = np.asarray(MFE, dtype=np.float64)
    if length.dtype.kind not in "iu":
        length, bad_len = _to_counts(np.asarray(length, dtype=np.float64))
    else:
        length = length.astype(np.int64, copy=False)
        bad_len = np.zeros(length.shape[0], dtype=bool)
    counts = []
    bad_count = bad_len
    for column in (NA, NC, NG, NT):
        column = np.asarray(column)
        if column.dtype.kind not in "iu":
            converted, bad = _to_counts(np.asarray(column, dtype=np.float64))
        else:
            converted = column.astype(np.int64, copy=False)
            bad = np.zeros(converted.shape[0], dtype=bool)
        counts.append(converted)
        bad_count = bad_count | bad
    na, nc, ng, nt = counts
    if length2 is None:
        length2_arr = None
        bad_l2 = np.zeros(length.shape[0], dtype=bool)
    else:
        length2_arr = np.asarray(length2)
        if length2_arr.dtype.kind not in "iu":
            length2_arr, bad_l2_exact = _to_counts(np.asarray(length2_arr, dtype=np.float64))
        else:
            length2_arr = length2_arr.astype(np.int64, copy=False)
            bad_l2_exact = np.zeros(length2_arr.shape[0], dtype=bool)
        bad_l2 = bad_l2_exact | (length2_arr != length)
    bad_L = length <= 0
    bad_sum = bad_count | ((na + nc + ng + nt) != length)
    bad_mfe = ~np.isfinite(mfe)
    bad = bad_L | bad_sum | bad_mfe | bad_l2
    reasons = {
        "L<=0": int(np.count_nonzero(bad_L)),
        "sum!=L": int(np.count_nonzero(bad_sum)),
        "mfe": int(np.count_nonzero(bad_mfe)),
        "length2": int(np.count_nonzero(bad_l2)),
    }
    return ~bad, reasons, length, mfe, na, nc, ng, nt


def _first_fields(path: str, separator: str) -> list[str]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        line = handle.readline()
    if line.startswith("\ufeff"):
        line = line[1:]
    line = line.rstrip("\r\n")
    return [field.strip() for field in line.split(separator)]


def _separator(path: str) -> str:
    suffix = Path(path).suffix.lower()
    if suffix == ".csv":
        return ","
    return "\t"


def _load_numeric(path: str, separator: str, skiprows: int, usecols: tuple[int, ...]) -> np.ndarray:
    try:
        data = np.loadtxt(
            path,
            delimiter=separator,
            usecols=usecols,
            skiprows=skiprows,
            dtype=np.float64,
            comments=None,
            ndmin=2,
            encoding="utf-8",
        )
    except ValueError as exc:
        raise ValueError(f"could not read numeric columns from {path}: {exc}") from exc
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def _arrays_from_columns(columns: dict[str, np.ndarray], mfe_float32: bool) -> dict:
    mfe = columns["MFE"]
    if mfe_float32:
        mfe = mfe.astype(np.float32)
    out = {
        "L": columns["L"],
        "MFE": mfe,
        "NA": columns["NA"],
        "NC": columns["NC"],
        "NG": columns["NG"],
        "NT": columns["NT"],
        "length2": columns.get("length2"),
    }
    return out


def _read_parquet(path: str, mfe_float32: bool) -> dict:
    try:
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise ValueError("reading parquet requires pyarrow") from exc
    names = list(pq.read_schema(path).names)
    index = {name: i for i, name in enumerate(names)}
    wanted = []
    for name in ("MFE", "NA", "NC", "NG", "NT"):
        if name not in index:
            raise ValueError(f"{path} is missing column {name}")
        wanted.append(name)
    if "L" in index:
        length_name = "L"
    elif "length" in index:
        length_name = "length"
    else:
        raise ValueError(f"{path} is missing column L (or length)")
    wanted.insert(0, length_name)
    if "length2" in index:
        wanted.append("length2")
    table = pq.read_table(path, columns=wanted)
    columns = {name: table.column(name).to_numpy(zero_copy_only=False) for name in wanted}
    columns["L"] = columns.pop(length_name)
    return _arrays_from_columns(columns, mfe_float32)


def read_table(path: str, mfe_float32: bool = False) -> dict:
    """Read MFE, L, NA, NC, NG, NT. ``length2`` is included when the file has it.

    ``.parquet`` uses pyarrow. ``.csv`` is comma-separated and ``.tsv`` / ``.txt``
    are tab-separated. A header is detected from the names above. A headerless
    file uses the production layout ``label length MFE length2 NA NC NG NT``
    and the label column is not loaded.
    """
    suffix = Path(path).suffix.lower()
    if suffix == ".parquet":
        return _read_parquet(path, mfe_float32)
    separator = _separator(path)
    fields = _first_fields(path, separator)
    if any(field in HEADER_TOKENS for field in fields):
        index = {name: i for i, name in enumerate(fields)}
        for name in ("MFE", "NA", "NC", "NG", "NT"):
            if name not in index:
                raise ValueError(f"{path} is missing column {name}")
        if "L" in index:
            length_at = index["L"]
        elif "length" in index:
            length_at = index["length"]
        else:
            raise ValueError(f"{path} is missing column L (or length)")
        usecols = [
            length_at,
            index["MFE"],
            index["NA"],
            index["NC"],
            index["NG"],
            index["NT"],
        ]
        names = ["L", "MFE", "NA", "NC", "NG", "NT"]
        if "length2" in index:
            usecols.append(index["length2"])
            names.append("length2")
        data = _load_numeric(path, separator, 1, tuple(usecols))
    else:
        names = ["L", "MFE", "length2", "NA", "NC", "NG", "NT"]
        data = _load_numeric(path, separator, 0, (1, 2, 3, 4, 5, 6, 7))
    columns = {name: data[:, i] for i, name in enumerate(names)}
    return _arrays_from_columns(columns, mfe_float32)


def _report_rows(n_read: int, n_dropped: int, reasons: dict, n_train: int, n_hold: int) -> None:
    print(
        f"read={n_read} dropped={n_dropped} "
        f"(L<=0: {reasons['L<=0']}, sum!=L: {reasons['sum!=L']}, "
        f"mfe: {reasons['mfe']}, length2: {reasons['length2']}) "
        f"train={n_train} holdout={n_hold}",
        file=sys.stderr,
        flush=True,
    )


def build_from_arrays(
    L, MFE, NA, NC, NG, NT, *, n_min, tau, max_depth, holdout_frac, seed, length2=None
):
    """Filter, hold out, and build. Returns ``(tree, rows_meta, X_train, y_train, X_hold, y_hold)``."""
    n_read = int(np.asarray(L).shape[0])
    keep, reasons, length, mfe, na, nc, ng, nt = filter_rows(
        L, MFE, NA, NC, NG, NT, length2=length2
    )
    del L, MFE, NA, NC, NG, NT, length2
    n_dropped = int(np.count_nonzero(~keep))
    length = length[keep]
    mfe = mfe[keep]
    na = na[keep]
    nc = nc[keep]
    ng = ng[keep]
    nt = nt[keep]
    n_keep = int(length.shape[0])
    if n_keep == 0:
        raise ValueError("no rows left after validation")
    train_idx, hold_idx = holdout_split(n_keep, holdout_frac, seed)
    if train_idx.size == 0:
        raise ValueError("training split is empty; lower holdout-frac")
    X_train = feat.features(length[train_idx], na[train_idx], nc[train_idx], ng[train_idx], nt[train_idx])
    X_hold = feat.features(length[hold_idx], na[hold_idx], nc[hold_idx], ng[hold_idx], nt[hold_idx])
    y_train = np.asarray(mfe[train_idx])
    y_hold = np.asarray(mfe[hold_idx])
    n_train = int(train_idx.size)
    n_hold = int(hold_idx.size)
    del length, mfe, na, nc, ng, nt, train_idx, hold_idx
    tree = build_tree(X_train, y_train, n_min, tau, max_depth)
    fill_holdout_stats(tree, X_hold, y_hold)
    rows = {
        "read": n_read,
        "dropped": n_dropped,
        "train": n_train,
        "holdout": n_hold,
        "dropped_reasons": reasons,
    }
    _report_rows(n_read, n_dropped, reasons, rows["train"], rows["holdout"])
    return tree, rows, X_train, y_train, X_hold, y_hold


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Build an adaptive MFE binning tree.")
    parser.add_argument("--input", required=True, help="CSV, TSV, or parquet table")
    parser.add_argument("--out-prefix", required=True, help="Writes PREFIX.tree.tsv and PREFIX.tree.meta.json")
    parser.add_argument("--n-min", type=int, default=2000)
    parser.add_argument("--tau", type=float, default=0.02)
    parser.add_argument("--max-depth", type=int, default=40)
    parser.add_argument("--holdout-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--mfe-float32", action="store_true", help="Store MFE as float32 before the float64 regression")
    args = parser.parse_args(argv)
    if args.n_min < 1:
        _die("--n-min must be >= 1")
    if args.tau < 0:
        _die("--tau must be >= 0")
    if args.max_depth < 0:
        _die("--max-depth must be >= 0")
    if not 0.0 <= args.holdout_frac <= 1.0:
        _die("--holdout-frac must be between 0 and 1")
    try:
        raw = read_table(args.input, mfe_float32=args.mfe_float32)
        tree, rows, X_train, y_train, X_hold, y_hold = build_from_arrays(
            raw.pop("L"),
            raw.pop("MFE"),
            raw.pop("NA"),
            raw.pop("NC"),
            raw.pop("NG"),
            raw.pop("NT"),
            length2=raw.pop("length2"),
            n_min=args.n_min,
            tau=args.tau,
            max_depth=args.max_depth,
            holdout_frac=args.holdout_frac,
            seed=args.seed,
        )
        del raw, X_train, y_train
    except ValueError as exc:
        _die(str(exc))
    z_stats = holdout_z_moments(tree, X_hold, y_hold)
    tsv_path = args.out_prefix + ".tree.tsv"
    meta_path = args.out_prefix + ".tree.meta.json"
    write_tree_tsv(tsv_path, tree)
    write_meta(
        meta_path,
        parameters={
            "n_min": int(args.n_min),
            "tau": float(args.tau),
            "max_depth": int(args.max_depth),
            "holdout_frac": float(args.holdout_frac),
            "seed": int(args.seed),
            "mfe_float32": bool(args.mfe_float32),
        },
        rows=rows,
        tree=tree,
    )
    print(format_summary(tree, z_stats))


if __name__ == "__main__":
    main()
