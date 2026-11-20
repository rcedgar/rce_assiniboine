#!/usr/bin/env python3
"""Shared feature map and vectorized lookup for MFE binning trees.

Build and assignment both call :func:`features`, so a threshold written at
build time is applied to the same float64 values later. The split convention
everywhere is ``x[axis] < threshold -> left`` and ``x >= threshold -> right``.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

# Natural log. x1 is the GC fraction, x2 the AT skew, x3 the GC skew.
FEATURE_DEFINITIONS: tuple[tuple[str, str], ...] = (
    ("x0", "log(L)"),
    ("x1", "(NC+NG)/L"),
    ("x2", "(NA-NT)/L"),
    ("x3", "(NG-NC)/L"),
)
N_FEATURES = 4
SPLIT_CONVENTION = "x[axis] < threshold -> left"

TREE_COLUMNS: tuple[str, ...] = (
    "node_id",
    "parent",
    "left",
    "right",
    "axis",
    "threshold",
    "depth",
    "is_leaf",
    "stop_reason",
    "n_train",
    "mean_train",
    "var_train",
    "skew_train",
    "exkurt_train",
    "q001_train",
    "q01_train",
    "q99_train",
    "q999_train",
    "n_hold",
    "mean_hold",
    "var_hold",
    "skew_hold",
    "exkurt_hold",
    "bbox_min_x0",
    "bbox_max_x0",
    "bbox_min_x1",
    "bbox_max_x1",
    "bbox_min_x2",
    "bbox_max_x2",
    "bbox_min_x3",
    "bbox_max_x3",
)


class Tree(NamedTuple):
    """Arrays indexed by node_id. Root is 0. Missing links and axes are -1.

    Leaf statistics and bounding boxes are NaN on internal nodes. ``n_hold``
    is 0 on a leaf that received no holdout rows, and NaN on internal nodes.
    ``threshold`` is NaN on leaves.
    """

    axis: np.ndarray
    thr: np.ndarray
    left: np.ndarray
    right: np.ndarray
    parent: np.ndarray
    is_leaf: np.ndarray
    depth: np.ndarray
    stop_reason: np.ndarray
    n_train: np.ndarray
    mean_train: np.ndarray
    var_train: np.ndarray
    skew_train: np.ndarray
    exkurt_train: np.ndarray
    q001_train: np.ndarray
    q01_train: np.ndarray
    q99_train: np.ndarray
    q999_train: np.ndarray
    n_hold: np.ndarray
    mean_hold: np.ndarray
    var_hold: np.ndarray
    skew_hold: np.ndarray
    exkurt_hold: np.ndarray
    bbox_min: np.ndarray
    bbox_max: np.ndarray


def features(L, NA, NC, NG, NT) -> np.ndarray:
    """Map sequence length and base counts to a float64 ``(n, 4)`` matrix.

    ``x0 = log(L)`` (natural log), ``x1 = (NC + NG) / L``,
    ``x2 = (NA - NT) / L``, ``x3 = (NG - NC) / L``.
    Inputs are cast to float64 before any arithmetic.
    """
    length = np.asarray(L, dtype=np.float64).reshape(-1)
    a = np.asarray(NA, dtype=np.float64).reshape(-1)
    c = np.asarray(NC, dtype=np.float64).reshape(-1)
    g = np.asarray(NG, dtype=np.float64).reshape(-1)
    t = np.asarray(NT, dtype=np.float64).reshape(-1)
    x0 = np.log(length)
    x1 = (c + g) / length
    x2 = (a - t) / length
    x3 = (g - c) / length
    return np.column_stack((x0, x1, x2, x3))


def _as_float_column(text: np.ndarray, index: int) -> np.ndarray:
    raw = text[:, index]
    cleaned = raw.copy()
    blank = (raw == "") | (raw == "NA") | (raw == "na")
    cleaned[blank] = "nan"
    return cleaned.astype(np.float64)


def load_tree(path: str) -> Tree:
    """Load a tree TSV written by the builder. Numeric fields come from numpy."""
    table = np.loadtxt(
        path,
        delimiter="\t",
        dtype=str,
        skiprows=1,
        comments=None,
        ndmin=2,
        encoding="utf-8",
    )
    if table.shape[1] != len(TREE_COLUMNS):
        raise ValueError(
            f"{path} has {table.shape[1]} columns, expected {len(TREE_COLUMNS)}"
        )
    node_id = table[:, 0].astype(np.int32)
    order = np.argsort(node_id, kind="mergesort")
    table = table[order]
    node_id = table[:, 0].astype(np.int32)
    if not np.array_equal(node_id, np.arange(node_id.shape[0], dtype=np.int32)):
        raise ValueError(f"{path} node_id values are not a contiguous range from 0")

    axis = table[:, 4].astype(np.int32)
    thr = _as_float_column(table, 5)
    left = table[:, 2].astype(np.int32)
    right = table[:, 3].astype(np.int32)
    parent = table[:, 1].astype(np.int32)
    is_leaf = table[:, 7].astype(np.int32).astype(bool)
    depth = table[:, 6].astype(np.int32)
    stop_reason = table[:, 8].astype(object)

    def fcol(index: int) -> np.ndarray:
        return _as_float_column(table, index)

    bbox_min = np.column_stack((fcol(23), fcol(25), fcol(27), fcol(29)))
    bbox_max = np.column_stack((fcol(24), fcol(26), fcol(28), fcol(30)))
    return Tree(
        axis=axis,
        thr=thr,
        left=left,
        right=right,
        parent=parent,
        is_leaf=is_leaf,
        depth=depth,
        stop_reason=stop_reason,
        n_train=fcol(9),
        mean_train=fcol(10),
        var_train=fcol(11),
        skew_train=fcol(12),
        exkurt_train=fcol(13),
        q001_train=fcol(14),
        q01_train=fcol(15),
        q99_train=fcol(16),
        q999_train=fcol(17),
        n_hold=fcol(18),
        mean_hold=fcol(19),
        var_hold=fcol(20),
        skew_hold=fcol(21),
        exkurt_hold=fcol(22),
        bbox_min=bbox_min,
        bbox_max=bbox_max,
    )


def points_in_bbox(tree: Tree, X: np.ndarray, leaf_ids: np.ndarray) -> np.ndarray:
    """True where each row of ``X`` lies inside its leaf's training bounding box."""
    X = np.asarray(X, dtype=np.float64)
    lo = tree.bbox_min[leaf_ids]
    hi = tree.bbox_max[leaf_ids]
    return np.all((X >= lo) & (X <= hi), axis=1)


def assign(tree: Tree, X: np.ndarray, *, return_in_bbox: bool = False):
    """Return the leaf node_id of each row of an ``(n, 4)`` feature matrix.

    The walk loops over tree depth only. A row moves to the right child when
    ``x[axis] >= threshold``, and to the left child otherwise.
    """
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2 or X.shape[1] != N_FEATURES:
        raise ValueError("X must have shape (n, 4)")
    n = X.shape[0]
    node = np.zeros(n, dtype=np.int32)
    if n == 0:
        if return_in_bbox:
            return node, np.zeros(0, dtype=bool)
        return node
    limit = int(np.max(tree.depth)) + 2
    for _ in range(limit):
        rows = np.flatnonzero(~tree.is_leaf[node])
        if rows.size == 0:
            break
        nid = node[rows]
        ax = tree.axis[nid]
        go_right = X[rows, ax] >= tree.thr[nid]
        node[rows] = tree.left[nid] + go_right.astype(np.int32)
    if not bool(np.all(tree.is_leaf[node])):
        raise RuntimeError("assignment did not reach leaves")
    if return_in_bbox:
        return node, points_in_bbox(tree, X, node)
    return node


def assign_counts(tree: Tree, L, NA, NC, NG, NT, *, return_in_bbox: bool = False):
    """Feature-map ``(L, NA, NC, NG, NT)`` and return leaf ids."""
    return assign(
        tree,
        features(L, NA, NC, NG, NT),
        return_in_bbox=return_in_bbox,
    )
