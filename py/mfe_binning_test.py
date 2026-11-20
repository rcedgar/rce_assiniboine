#!/usr/bin/env python3
"""Synthetic tests for the adaptive MFE partition."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import mfe_binning_build_bins as build
import mfe_binning_features as feat

N_MIN = 400
TAU = 0.02
MAX_DEPTH = 40


def make_smooth(n: int = 200_000, seed: int = 1):
    """Log-uniform lengths, random compositions, MFE linear in L and L*GC."""
    rng = np.random.default_rng(seed)
    log_L = rng.uniform(np.log(20.0), np.log(2000.0), n)
    L = np.clip(np.rint(np.exp(log_L)), 20, 2000).astype(np.int64)
    # Sequential draws so NA+NC+NG+NT == L exactly, with no per-base remainder fixup.
    NA = rng.integers(0, L + 1)
    NC = rng.integers(0, L - NA + 1)
    NG = rng.integers(0, L - NA - NC + 1)
    NT = L - NA - NC - NG
    gc = (NC + NG) / L.astype(np.float64)
    mfe = -0.3 * L + -0.05 * L * gc + rng.normal(0.0, 25.0, n)
    return L, NA, NC, NG, NT, mfe


def make_ties(seed: int = 2):
    """Few distinct (L, composition) values, repeated."""
    rng = np.random.default_rng(seed)
    lengths = (40, 80, 120)
    bases = (
        (1, 1, 1, 1),
        (2, 1, 1, 0),
        (3, 1, 0, 0),
        (1, 2, 1, 0),
        (0, 1, 1, 2),
    )
    rows_L = []
    rows_counts = []
    for length in lengths:
        for base in bases:
            scale = length // sum(base)
            counts = tuple(part * scale for part in base)
            rows_L.append(np.full(200, length, dtype=np.int64))
            rows_counts.append(np.tile(np.array(counts, dtype=np.int64), (200, 1)))
    L = np.concatenate(rows_L)
    counts = np.vstack(rows_counts)
    mfe = -0.3 * L + rng.normal(0.0, 1.0, size=L.shape[0])
    return L, counts[:, 0], counts[:, 1], counts[:, 2], counts[:, 3], mfe


def assert_tree_structure(tree: feat.Tree, n_min: int) -> None:
    n_nodes = int(tree.is_leaf.shape[0])
    assert int(np.sum(tree.parent == -1)) == 1
    assert int(tree.parent[0]) == -1
    n_leaves = int(np.count_nonzero(tree.is_leaf))
    assert n_leaves == (n_nodes - n_leaves) + 1
    seen = np.zeros(n_nodes, dtype=bool)
    stack = [0]
    while stack:
        node_id = stack.pop()
        assert not seen[node_id]
        seen[node_id] = True
        if tree.is_leaf[node_id]:
            assert int(tree.left[node_id]) == -1
            assert int(tree.right[node_id]) == -1
            assert int(tree.axis[node_id]) == -1
            assert tree.stop_reason[node_id] in {"min_n", "max_depth", "bias_ok", "no_split"}
            assert int(tree.n_train[node_id]) >= n_min
        else:
            left = int(tree.left[node_id])
            right = int(tree.right[node_id])
            assert right == left + 1
            assert int(tree.parent[left]) == node_id
            assert int(tree.parent[right]) == node_id
            assert int(tree.axis[node_id]) >= 0
            assert np.isfinite(tree.thr[node_id])
            stack.append(right)
            stack.append(left)
    assert bool(seen.all())
    sizes = np.zeros(n_nodes, dtype=np.int64)
    sizes[tree.is_leaf] = np.rint(tree.n_train[tree.is_leaf]).astype(np.int64)
    for node_id in range(n_nodes - 1, -1, -1):
        if tree.is_leaf[node_id]:
            continue
        left = int(tree.left[node_id])
        right = int(tree.right[node_id])
        sizes[node_id] = sizes[left] + sizes[right]
        assert sizes[left] >= n_min
        assert sizes[right] >= n_min


@pytest.fixture(scope="module")
def smooth():
    L, NA, NC, NG, NT, mfe = make_smooth()
    assert np.all(NA + NC + NG + NT == L)
    tree, rows, X_train, y_train, X_hold, y_hold = build.build_from_arrays(
        L,
        mfe,
        NA,
        NC,
        NG,
        NT,
        n_min=N_MIN,
        tau=TAU,
        max_depth=MAX_DEPTH,
        holdout_frac=0.2,
        seed=1,
    )
    return {
        "tree": tree,
        "rows": rows,
        "X_train": X_train,
        "y_train": y_train,
        "X_hold": X_hold,
        "y_hold": y_hold,
        "L": L,
        "NA": NA,
        "NC": NC,
        "NG": NG,
        "NT": NT,
    }


def test_features_and_synthetic_tree(smooth):
    sample = feat.features([10], [1], [2], [3], [4])
    assert sample.shape == (1, 4)
    assert sample.dtype == np.float64
    assert sample[0, 0] == np.log(10.0)
    assert sample[0, 1] == 0.5
    assert sample[0, 2] == (1.0 - 4.0) / 10.0
    assert sample[0, 3] == (3.0 - 2.0) / 10.0

    tree = smooth["tree"]
    rows = smooth["rows"]
    assert rows["read"] == rows["dropped"] + rows["train"] + rows["holdout"]
    assert rows["dropped"] == 0
    assert rows["train"] == smooth["X_train"].shape[0]
    assert int(np.count_nonzero(~tree.is_leaf)) >= 1
    reasons = tree.stop_reason[tree.is_leaf]
    counts = {name: int(np.sum(reasons == name)) for name in ("bias_ok", "min_n", "max_depth", "no_split")}
    assert counts["bias_ok"] >= 1, counts


def test_training_rows_match_leaf_counts(smooth):
    tree = smooth["tree"]
    leaf_ids = feat.assign(tree, smooth["X_train"])
    assert bool(np.all(tree.is_leaf[leaf_ids]))
    counts = np.bincount(leaf_ids, minlength=tree.is_leaf.shape[0])
    assert np.all(counts[~tree.is_leaf] == 0)
    expected = np.rint(tree.n_train[tree.is_leaf]).astype(np.int64)
    assert np.array_equal(counts[tree.is_leaf], expected)
    ids_again, inside = feat.assign(tree, smooth["X_train"], return_in_bbox=True)
    assert np.array_equal(ids_again, leaf_ids)
    assert bool(np.all(inside))


def test_leaf_min_size_and_structure(smooth):
    assert_tree_structure(smooth["tree"], N_MIN)


def test_tree_tsv_round_trip(smooth, tmp_path):
    tree = smooth["tree"]
    path = tmp_path / "round.tree.tsv"
    build.write_tree_tsv(str(path), tree)
    loaded = feat.load_tree(str(path))
    original = feat.assign(tree, smooth["X_train"])
    restored = feat.assign(loaded, smooth["X_train"])
    assert np.array_equal(original, restored)
    internal = ~tree.is_leaf
    assert np.array_equal(tree.thr[internal], loaded.thr[internal])
    assert np.array_equal(tree.left, loaded.left)
    assert np.array_equal(tree.axis, loaded.axis)


def test_tie_values_have_no_empty_children():
    L, NA, NC, NG, NT, mfe = make_ties()
    assert np.all(NA + NC + NG + NT == L)
    n_min = 50
    tree, rows, X_train, _y_train, _X_hold, _y_hold = build.build_from_arrays(
        L,
        mfe,
        NA,
        NC,
        NG,
        NT,
        n_min=n_min,
        tau=TAU,
        max_depth=MAX_DEPTH,
        holdout_frac=0.2,
        seed=2,
    )
    assert rows["train"] == X_train.shape[0]
    assert int(np.count_nonzero(~tree.is_leaf)) >= 1
    assert_tree_structure(tree, n_min)
    # Distinct feature rows are far fewer than the training set.
    unique = np.unique(np.round(X_train, decimals=12), axis=0)
    assert unique.shape[0] < 30


def test_bias_ok_leaves_meet_tau(smooth):
    tree = smooth["tree"]
    leaf_ids = feat.assign(tree, smooth["X_train"])
    bias_leaves = np.flatnonzero((tree.stop_reason == "bias_ok") & tree.is_leaf)
    assert bias_leaves.size >= 1
    for leaf_id in bias_leaves:
        rows = np.flatnonzero(leaf_ids == int(leaf_id)).astype(np.int32)
        decision = build.evaluate_node(
            smooth["X_train"],
            smooth["y_train"],
            rows,
            depth=0,
            n_min=N_MIN,
            tau=TAU,
            max_depth=MAX_DEPTH,
        )
        assert decision.stop_reason == "bias_ok"
        assert decision.s2 > 0.0
        assert decision.bias_var / decision.s2 < TAU


def test_stable_threshold_round_trip():
    lo = np.float64(1.0)
    hi = np.nextafter(lo, np.float64(np.inf))
    thr = build.stable_threshold(lo, hi)
    parsed = np.float64(format(float(thr), ".17g"))
    assert lo < thr <= hi
    assert parsed == thr
    assert lo < hi
    assert lo < thr
    assert hi >= thr

    midpoint = build.stable_threshold(np.float64(1.0), np.float64(3.0))
    assert midpoint == np.float64(2.0)
    assert np.float64(format(float(midpoint), ".17g")) == midpoint
    assert build.find_cut(np.ones(100, dtype=np.float64), n_min=10) is None


def test_headerless_production_layout_and_named_csv(tmp_path):
    headerless = tmp_path / "merged_compos_mfe.tsv"
    headerless.write_text(
        "\n".join(
            [
                "lab\t10\t-1.5\t10\t1\t2\t3\t4",
                "bad\t10\t-1.0\t10\t1\t2\t3\t5",
                "nanm\t8\tnan\t8\t2\t2\t2\t2",
                "l2\t8\t-2\t7\t2\t2\t2\t2",
                "ok\t8\t-3.25\t8\t2\t2\t2\t2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    raw = build.read_table(str(headerless))
    keep, reasons, length, mfe, na, nc, ng, nt = build.filter_rows(
        raw["L"], raw["MFE"], raw["NA"], raw["NC"], raw["NG"], raw["NT"], length2=raw["length2"]
    )
    assert reasons == {"L<=0": 0, "sum!=L": 1, "mfe": 1, "length2": 1}
    assert int(np.count_nonzero(keep)) == 2
    assert np.array_equal(length[keep], np.array([10, 8]))
    assert mfe[keep][1] == pytest.approx(-3.25)

    named = tmp_path / "named.csv"
    named.write_text("label,length,MFE,length2,NA,NC,NG,NT\nok,12,-5,12,3,3,3,3\n", encoding="utf-8")
    named_raw = build.read_table(str(named))
    assert named_raw["length2"] is not None
    assert int(named_raw["L"][0]) == 12
    assert float(named_raw["MFE"][0]) == -5.0
    x = feat.features(named_raw["L"], named_raw["NA"], named_raw["NC"], named_raw["NG"], named_raw["NT"])
    assert x.shape == (1, 4)


def test_cli_writes_tree_and_meta(tmp_path):
    src = tmp_path / "small.tsv"
    L, NA, NC, NG, NT, mfe = make_smooth(n=2000, seed=3)
    with src.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("MFE\tL\tNA\tNC\tNG\tNT\n")
        for i in range(L.shape[0]):
            handle.write(
                f"{mfe[i]:.8g}\t{int(L[i])}\t{int(NA[i])}\t{int(NC[i])}\t{int(NG[i])}\t{int(NT[i])}\n"
            )
    prefix = str(tmp_path / "run1")
    build.main(
        [
            "--input",
            str(src),
            "--out-prefix",
            prefix,
            "--n-min",
            "50",
            "--tau",
            "0.02",
            "--holdout-frac",
            "0.2",
            "--seed",
            "1",
        ]
    )
    loaded = feat.load_tree(prefix + ".tree.tsv")
    meta = (tmp_path / "run1.tree.meta.json").read_text(encoding="utf-8")
    assert "x[axis] < threshold -> left" in meta
    assert int(np.count_nonzero(loaded.is_leaf)) >= 1
    assert loaded.parent[0] == -1
