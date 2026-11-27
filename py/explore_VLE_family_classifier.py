#!/usr/bin/env python3
"""Find circle features that separate VLE families from each other.

DECOY rows are discarded. Each remaining row is one VLE circle, and the
family is the second field of the label. The question is which features
support a partition of the families, not a VLE-versus-DECOY call.

Epsilon-squared is the tie-corrected Kruskal-Wallis H divided by n - 1.
Oriented AUC is max(AUC, 1 - AUC). A subset tree keeps a split when its
best single-feature oriented AUC is at least 0.70. The ridge model is
only a ceiling for a later classifier.
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from typing import NamedTuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_INPUT = _ROOT / "big_output" / "controls_circle_features_filtered.tsv"
_DEFAULT_OUTPUT = _ROOT / "explore_VLE_family_classifier"
_LENGTH = "nr_bases"
_NEAR_CONSTANT_FRACTION = 0.99
_UNSTABLE_N = 20
_SPLIT_MIN = 0.70
_CLEAR_AUC = 0.80
_TOP_N = 8
_BOX_N = 6
_LENGTH_RHO_MAX = 0.30
_PCA_PER_FAMILY = 300
_FOLDS = 5
_RIDGE = 1.0
_SEED = 1
_COLOR_RAW = "#d55e00"
_COLOR_ADJ = "#0072b2"
_COLORS = (
    "#0072b2",
    "#d55e00",
    "#009e73",
    "#cc79a7",
    "#e69f00",
    "#56b4e9",
    "#f0e442",
    "#000000",
)
_IDENTITY_NAMES = {
    "bulge_to_stem_ratio",
    "bKeto",
    "uKeto",
    "pS",
    "pH",
    "pI",
    "pB",
    "pF",
    "pRod",
    "pLongestRod",
}


class _FeatureScore(NamedTuple):
    name: str
    epsilon2: float
    h: float
    max_ovr: float
    extreme_family: str
    extreme_direction: str
    epsilon2_adj: float | None
    h_adj: float | None
    max_ovr_adj: float | None
    extreme_family_adj: str | None
    extreme_direction_adj: str | None


class _Pair(NamedTuple):
    family_a: str
    family_b: str
    n_a: int
    n_b: int
    unstable: bool
    feature: str
    oriented_auc: float
    family_a_direction: str
    nr_bases_oriented_auc: float
    feature_adj: str
    oriented_auc_adj: float
    family_a_direction_adj: str


class _Ovr(NamedTuple):
    family: str
    n: int
    unstable: bool
    feature: str
    oriented_auc: float
    direction: str
    feature_adj: str
    oriented_auc_adj: float
    direction_adj: str


class _Split(NamedTuple):
    view: str
    depth: int
    parent: str
    feature: str
    oriented_auc: float
    threshold: float
    lower_families: tuple[str, ...]
    higher_families: tuple[str, ...]
    n_lower: int
    n_higher: int
    median_lower: float
    median_higher: float
    unstable: bool
    expanded: str


class _Ridge(NamedTuple):
    name: str
    balanced_accuracy: float
    recalls: list[float]
    confusion: np.ndarray


class _Pick(NamedTuple):
    oriented: float
    epsilon2: float
    name: str
    direction: str
    column: int


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _expect(condition: bool, message: str) -> None:
    if not condition:
        _die("ERROR: self-check failed: " + message)


def _fmt(value: float | None) -> str:
    if value is None:
        return "NA"
    if not math.isfinite(value):
        return "inf" if value > 0.0 else "-inf"
    return "%.6g" % value


def _yes(flag: bool) -> str:
    return "yes" if flag else "no"


def _join(names: list[str] | tuple[str, ...]) -> str:
    return ",".join(names)


def _is_header(first_field: str) -> bool:
    return first_field in ("#label", "label")


def _transform_of(name: str) -> str:
    if name.startswith("pct_") or "_pct_" in name or name in _IDENTITY_NAMES:
        return "identity"
    return "log1p"


def _family_of(label: str, line_number: int) -> str | None:
    parts = label.split("_", 2)
    if len(parts) != 3 or parts[0] not in ("VLE", "DECOY") or not parts[1] or not parts[2]:
        _die(
            "ERROR: line %d label '%s' is not VLE_family_... or DECOY_family_..."
            % (line_number, label)
        )
    if parts[0] == "DECOY":
        return None
    return parts[1]


def _canonical_feature_names() -> list[str]:
    py_dir = str(Path(__file__).resolve().parent)
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)
    from dotbracket_circular_features import extract_circular_rna_features

    names = sorted(extract_circular_rna_features("").keys())
    if not names:
        _die("ERROR: extract_circular_rna_features returned no features")
    return names


def _load(path: Path) -> tuple[list[str], np.ndarray, list[str], int]:
    """Return feature names, VLE values, the family of each VLE row, and the DECOY count."""
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            first = handle.readline()
            if first == "":
                _die("ERROR: %s is empty" % path)
            first_fields = first.lstrip("\ufeff").rstrip("\r\n").split("\t")
            if _is_header(first_fields[0]):
                feature_names = first_fields[1:]
                pending_first = None
            else:
                feature_names = _canonical_feature_names()
                pending_first = first_fields
            if not feature_names:
                _die("ERROR: %s has no feature columns" % path)
            if len(feature_names) != len(set(feature_names)):
                _die("ERROR: %s repeats a feature name" % path)
            if _LENGTH not in feature_names:
                _die("ERROR: %s has no %s column" % (path, _LENGTH))
            width = 1 + len(feature_names)
            columns: list[list[float]] = [[] for _name in feature_names]
            families: list[str] = []
            n_decoy = 0

            def consume(line_number: int, fields: list[str]) -> None:
                nonlocal n_decoy
                if len(fields) != width:
                    _die(
                        "ERROR: %s:%d has %d columns, expected %d"
                        % (path, line_number, len(fields), width)
                    )
                family = _family_of(fields[0], line_number)
                if family is None:
                    n_decoy += 1
                    return
                families.append(family)
                for index, raw in enumerate(fields[1:]):
                    try:
                        value = float(raw)
                    except ValueError:
                        _die(
                            "ERROR: %s:%d %s is not a number"
                            % (path, line_number, feature_names[index])
                        )
                    if not math.isfinite(value):
                        _die(
                            "ERROR: %s:%d %s is not finite"
                            % (path, line_number, feature_names[index])
                        )
                    columns[index].append(value)

            if pending_first is not None:
                consume(1, pending_first)
            next_report = 100000
            seen = 1 if pending_first is not None else 1
            for line_number, line in enumerate(handle, start=2):
                text = line.rstrip("\r\n")
                if text == "":
                    continue
                consume(line_number, text.split("\t"))
                seen += 1
                if seen >= next_report:
                    print("scanned %d rows" % seen, file=sys.stderr, flush=True)
                    next_report += 100000
    except OSError as exc:
        _die("ERROR: cannot read '%s': %s" % (path, exc.strerror or exc))
    if not families:
        _die("ERROR: no VLE rows in %s" % path)
    values = np.column_stack([np.asarray(column, dtype=np.float64) for column in columns])
    return feature_names, values, families, n_decoy


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """1-based ranks. Tied values share the mean of the ranks they span."""
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(values.size, dtype=np.float64)
    start = 0
    n = values.size
    while start < n:
        end = start + 1
        while end < n and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[order[start:end]] = 0.5 * ((start + 1) + end)
        start = end
    return ranks


def _auc_from_ranks(ranks: np.ndarray, is_a: np.ndarray) -> float:
    """P(A > B) + 0.5 P(tie). is_a marks group A."""
    n_a = int(is_a.sum())
    n_b = int(ranks.size - n_a)
    if n_a == 0 or n_b == 0:
        _die("ERROR: AUC group is empty")
    rank_sum = float(ranks[is_a].sum())
    return (rank_sum - n_a * (n_a + 1) / 2.0) / float(n_a * n_b)


def _orient(auc: float) -> tuple[float, str]:
    """Return oriented AUC and whether group A tends to be larger."""
    if auc > 0.5:
        return auc, "higher"
    if auc < 0.5:
        return 1.0 - auc, "lower"
    return 0.5, "tie"


def _tie_factor(values: np.ndarray) -> float:
    ordered = np.sort(values, kind="mergesort")
    n = int(ordered.size)
    acc = 0.0
    start = 0
    while start < n:
        end = start + 1
        while end < n and ordered[end] == ordered[start]:
            end += 1
        t = end - start
        if t > 1:
            acc += float(t * t * t - t)
        start = end
    denom = float(n * n * n - n)
    if denom <= 0.0 or acc <= 0.0:
        return 1.0
    factor = 1.0 - acc / denom
    if factor <= 0.0:
        return 0.0
    return factor


def _kruskal(values: np.ndarray, groups: np.ndarray, n_groups: int) -> tuple[float, float]:
    """Tie-corrected Kruskal-Wallis H and epsilon-squared H / (n - 1)."""
    n = int(values.size)
    if n < 2 or np.all(values == values[0]):
        return 0.0, 0.0
    ranks = _average_ranks(values)
    h_sum = 0.0
    for group in range(n_groups):
        mask = groups == group
        n_group = int(mask.sum())
        if n_group == 0:
            _die("ERROR: Kruskal-Wallis group is empty")
        rank_sum = float(ranks[mask].sum())
        h_sum += (rank_sum * rank_sum) / float(n_group)
    h = (12.0 / (n * (n + 1.0))) * h_sum - 3.0 * (n + 1.0)
    factor = _tie_factor(values)
    if factor == 0.0:
        return 0.0, 0.0
    h = h / factor
    if not math.isfinite(h):
        _die("ERROR: Kruskal-Wallis statistic is not finite")
    if h < 0.0:
        if h > -1e-6:
            h = 0.0
        else:
            _die("ERROR: Kruskal-Wallis statistic is negative (%s)" % _fmt(h))
    epsilon2 = h / float(n - 1)
    if epsilon2 > 1.0 + 1e-6:
        _die("ERROR: epsilon-squared is above 1 (%s)" % _fmt(epsilon2))
    if epsilon2 > 1.0:
        epsilon2 = 1.0
    return h, epsilon2


def _residualize(values: np.ndarray, log_length: np.ndarray) -> np.ndarray:
    design = np.column_stack((np.ones(values.shape[0]), log_length))
    coefficient, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    if not np.all(np.isfinite(coefficient)):
        _die("ERROR: length residual is not finite")
    return values - design @ coefficient


def _residualize_fold(
    train: np.ndarray,
    test: np.ndarray,
    log_train: np.ndarray,
    log_test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    design = np.column_stack((np.ones(train.shape[0]), log_train))
    coefficient, _, _, _ = np.linalg.lstsq(design, train, rcond=None)
    if not np.all(np.isfinite(coefficient)):
        _die("ERROR: fold length residual is not finite")
    test_design = np.column_stack((np.ones(test.shape[0]), log_test))
    return train - design @ coefficient, test - test_design @ coefficient


def _modal_fraction(values: np.ndarray) -> tuple[float, float]:
    _unique, counts = np.unique(values, return_counts=True)
    mode_index = int(np.argmax(counts))
    return float(_unique[mode_index]), float(counts[mode_index] / values.size)


def _near_constant(
    names: list[str], values: np.ndarray
) -> list[tuple[str, float, float]]:
    found: list[tuple[str, float, float]] = []
    for index, name in enumerate(names):
        mode, fraction = _modal_fraction(values[:, index])
        if fraction >= _NEAR_CONSTANT_FRACTION:
            found.append((name, mode, fraction))
    found.sort(key=lambda item: (-item[2], item[0]))
    return found


def _duplicate_groups(values: np.ndarray, names: list[str]) -> list[tuple[str, ...]]:
    used: set[int] = set()
    groups: list[tuple[str, ...]] = []
    for i, name in enumerate(names):
        if i in used:
            continue
        group = [name]
        for j in range(i + 1, len(names)):
            if j in used:
                continue
            if np.array_equal(values[:, i], values[:, j]):
                group.append(names[j])
                used.add(j)
        if len(group) > 1:
            groups.append(tuple(sorted(group)))
            used.add(i)
    groups.sort()
    return groups


def _bipartitions(codes: list[int]) -> list[tuple[list[int], list[int]]]:
    """Every nontrivial split. codes[0] stays on side A so complements are not repeated."""
    if len(codes) < 2:
        return []
    others = codes[1:]
    found: list[tuple[list[int], list[int]]] = []
    for mask in range(1 << len(others)):
        side_a = [codes[0]]
        side_b: list[int] = []
        for bit, code in enumerate(others):
            if mask & (1 << bit):
                side_a.append(code)
            else:
                side_b.append(code)
        if not side_b:
            continue
        found.append((side_a, side_b))
    return found


def _leaf_order(dist: np.ndarray) -> list[int]:
    """Average-linkage leaf order. The child with the smaller leaf index goes left."""
    k = int(dist.shape[0])
    if k == 0:
        return []
    clusters: list[list[int]] = [[i] for i in range(k)]
    current = dist.astype(np.float64).copy()
    np.fill_diagonal(current, np.inf)
    while len(clusters) > 1:
        flat = int(np.argmin(current))
        i, j = divmod(flat, current.shape[0])
        if i > j:
            i, j = j, i
        if not math.isfinite(current[i, j]):
            _die("ERROR: clustering distance is not finite")
        n_i = len(clusters[i])
        n_j = len(clusters[j])
        updated: list[float] = []
        for other in range(len(clusters)):
            if other == i or other == j:
                continue
            updated.append((n_i * current[i, other] + n_j * current[j, other]) / float(n_i + n_j))
        left = clusters[i]
        right = clusters[j]
        if min(left) > min(right):
            left, right = right, left
        merged = left + right
        keep = [other for other in range(len(clusters)) if other != i and other != j]
        clusters = [clusters[other] for other in keep]
        clusters.append(merged)
        size = len(clusters)
        rebuilt = np.full((size, size), np.inf, dtype=np.float64)
        for a in range(size - 1):
            for b in range(a + 1, size - 1):
                distance = current[keep[a], keep[b]]
                rebuilt[a, b] = distance
                rebuilt[b, a] = distance
        for a, distance in enumerate(updated):
            rebuilt[a, size - 1] = distance
            rebuilt[size - 1, a] = distance
        current = rebuilt
    return clusters[0]


def _better(pick: _Pick, best: _Pick | None) -> bool:
    if best is None:
        return True
    return (pick.oriented, pick.epsilon2, best.name) > (
        best.oriented,
        best.epsilon2,
        pick.name,
    )


def _best_feature(
    ranks_by_column: list[np.ndarray],
    columns: list[int],
    is_a: np.ndarray,
    epsilon2: np.ndarray,
    names: list[str],
) -> _Pick:
    best: _Pick | None = None
    for slot, column in enumerate(columns):
        auc = _auc_from_ranks(ranks_by_column[slot], is_a)
        oriented, direction = _orient(auc)
        pick = _Pick(oriented, float(epsilon2[column]), names[column], direction, column)
        if _better(pick, best):
            best = pick
    if best is None:
        _die("ERROR: no feature to score")
    return best


def _score_features(
    names: list[str],
    values: np.ndarray,
    residual: np.ndarray,
    groups: np.ndarray,
    family_names: list[str],
    length_index: int,
    eligible: list[int],
    eligible_adj: list[int],
    near_names: set[str],
) -> tuple[list[_FeatureScore], list[_Ovr]]:
    n_groups = len(family_names)
    scores: list[_FeatureScore] = []
    ovr_best: list[tuple[float, int, str, str]] = [
        (-1.0, -1, "", "") for _group in range(n_groups)
    ]
    ovr_best_adj: list[tuple[float, int, str, str]] = [
        (-1.0, -1, "", "") for _group in range(n_groups)
    ]
    ovr_raw: dict[int, list[tuple[float, str]]] = {column: [] for column in eligible}
    ovr_adj: dict[int, list[tuple[float, str]]] = {column: [] for column in eligible_adj}
    h_raw = np.zeros(len(names), dtype=np.float64)
    eps_raw = np.zeros(len(names), dtype=np.float64)
    h_adj = np.full(len(names), np.nan, dtype=np.float64)
    eps_adj = np.full(len(names), np.nan, dtype=np.float64)

    for column, name in enumerate(names):
        h, epsilon2 = _kruskal(values[:, column], groups, n_groups)
        h_raw[column] = h
        eps_raw[column] = epsilon2
        if column == length_index or name in near_names:
            continue
        h, epsilon2 = _kruskal(residual[:, column], groups, n_groups)
        h_adj[column] = h
        eps_adj[column] = epsilon2

    for column in eligible:
        ranks = _average_ranks(values[:, column])
        found: list[tuple[float, str]] = []
        for group, family in enumerate(family_names):
            oriented, direction = _orient(_auc_from_ranks(ranks, groups == group))
            found.append((oriented, direction))
            current = ovr_best[group]
            if (oriented, eps_raw[column], current[2]) > (
                current[0],
                eps_raw[current[1]] if current[1] >= 0 else -1.0,
                name,
            ):
                ovr_best[group] = (oriented, column, name, direction)
        ovr_raw[column] = found

    for column in eligible_adj:
        ranks = _average_ranks(residual[:, column])
        found = []
        for group, family in enumerate(family_names):
            oriented, direction = _orient(_auc_from_ranks(ranks, groups == group))
            found.append((oriented, direction))
            current = ovr_best_adj[group]
            if (oriented, eps_adj[column], current[2]) > (
                current[0],
                eps_adj[current[1]] if current[1] >= 0 else -1.0,
                names[column],
            ):
                ovr_best_adj[group] = (oriented, column, names[column], direction)
        ovr_adj[column] = found

    for column, name in enumerate(names):
        extreme_auc = -1.0
        extreme_family = family_names[0]
        extreme_direction = "tie"
        if column in ovr_raw:
            for group, (oriented, direction) in enumerate(ovr_raw[column]):
                if oriented > extreme_auc or (
                    oriented == extreme_auc and family_names[group] < extreme_family
                ):
                    extreme_auc = oriented
                    extreme_family = family_names[group]
                    extreme_direction = direction
        if column not in ovr_raw:
            extreme_auc = 0.5
            extreme_family = "NA"
            extreme_direction = "tie"
        if column == length_index or name in near_names:
            scores.append(
                _FeatureScore(
                    name=name,
                    epsilon2=float(eps_raw[column]),
                    h=float(h_raw[column]),
                    max_ovr=extreme_auc,
                    extreme_family=extreme_family,
                    extreme_direction=extreme_direction,
                    epsilon2_adj=None,
                    h_adj=None,
                    max_ovr_adj=None,
                    extreme_family_adj=None,
                    extreme_direction_adj=None,
                )
            )
            continue
        extreme_auc_adj = -1.0
        extreme_family_adj = family_names[0]
        extreme_direction_adj = "tie"
        if column in ovr_adj:
            for group, (oriented, direction) in enumerate(ovr_adj[column]):
                if oriented > extreme_auc_adj or (
                    oriented == extreme_auc_adj and family_names[group] < extreme_family_adj
                ):
                    extreme_auc_adj = oriented
                    extreme_family_adj = family_names[group]
                    extreme_direction_adj = direction
        else:
            extreme_auc_adj = 0.5
            extreme_family_adj = "NA"
            extreme_direction_adj = "tie"
        scores.append(
            _FeatureScore(
                name=name,
                epsilon2=float(eps_raw[column]),
                h=float(h_raw[column]),
                max_ovr=extreme_auc,
                extreme_family=extreme_family,
                extreme_direction=extreme_direction,
                epsilon2_adj=float(eps_adj[column]),
                h_adj=float(h_adj[column]),
                max_ovr_adj=extreme_auc_adj,
                extreme_family_adj=extreme_family_adj,
                extreme_direction_adj=extreme_direction_adj,
            )
        )
    scores.sort(key=lambda score: (-score.epsilon2, score.name))

    one_vs_rest: list[_Ovr] = []
    counts = np.bincount(groups, minlength=n_groups)
    for group, family in enumerate(family_names):
        raw = ovr_best[group]
        adj = ovr_best_adj[group]
        if raw[1] < 0 or adj[1] < 0:
            _die("ERROR: one-versus-rest search found no feature")
        one_vs_rest.append(
            _Ovr(
                family=family,
                n=int(counts[group]),
                unstable=int(counts[group]) < _UNSTABLE_N,
                feature=raw[2],
                oriented_auc=float(raw[0]),
                direction=raw[3],
                feature_adj=adj[2],
                oriented_auc_adj=float(adj[0]),
                direction_adj=adj[3],
            )
        )
    one_vs_rest.sort(key=lambda row: (-row.oriented_auc, row.family))
    return scores, one_vs_rest


def _pairwise(
    names: list[str],
    values: np.ndarray,
    residual: np.ndarray,
    groups: np.ndarray,
    family_names: list[str],
    counts: np.ndarray,
    length_index: int,
    eligible: list[int],
    eligible_adj: list[int],
    epsilon2: np.ndarray,
    epsilon2_adj: np.ndarray,
) -> list[_Pair]:
    pairs: list[_Pair] = []
    n_groups = len(family_names)
    for left in range(n_groups):
        for right in range(left + 1, n_groups):
            mask = (groups == left) | (groups == right)
            is_a = groups[mask] == left
            raw_ranks = [_average_ranks(values[mask, column]) for column in eligible]
            raw = _best_feature(raw_ranks, eligible, is_a, epsilon2, names)
            length_ranks = _average_ranks(values[mask, length_index])
            length_auc, _direction = _orient(_auc_from_ranks(length_ranks, is_a))
            adj_ranks = [_average_ranks(residual[mask, column]) for column in eligible_adj]
            adj = _best_feature(adj_ranks, eligible_adj, is_a, epsilon2_adj, names)
            pairs.append(
                _Pair(
                    family_a=family_names[left],
                    family_b=family_names[right],
                    n_a=int(counts[left]),
                    n_b=int(counts[right]),
                    unstable=int(counts[left]) < _UNSTABLE_N or int(counts[right]) < _UNSTABLE_N,
                    feature=raw.name,
                    oriented_auc=raw.oriented,
                    family_a_direction=raw.direction,
                    nr_bases_oriented_auc=length_auc,
                    feature_adj=adj.name,
                    oriented_auc_adj=adj.oriented,
                    family_a_direction_adj=adj.direction,
                )
            )
    pairs.sort(key=lambda row: (-row.oriented_auc, row.family_a, row.family_b))
    return pairs


def _best_split(
    codes: list[int],
    view: str,
    depth: int,
    names: list[str],
    matrix: np.ndarray,
    groups: np.ndarray,
    family_names: list[str],
    counts: np.ndarray,
    columns: list[int],
    epsilon2: np.ndarray,
) -> tuple[_Split, list[int], list[int]]:
    present = np.isin(groups, np.asarray(codes, dtype=np.int32))
    sub_groups = groups[present]
    sub_values = matrix[present]
    ranks = [_average_ranks(sub_values[:, column]) for column in columns]
    best: _Pick | None = None
    best_sides: tuple[list[int], list[int]] | None = None
    for side_a, side_b in _bipartitions(codes):
        is_a = np.isin(sub_groups, np.asarray(side_a, dtype=np.int32))
        pick = _best_feature(ranks, columns, is_a, epsilon2, names)
        if _better(pick, best):
            best = pick
            best_sides = (list(side_a), list(side_b))
    if best is None or best_sides is None:
        _die("ERROR: subset search found no split")
    side_a, side_b = best_sides
    if best.direction == "higher":
        higher_codes, lower_codes = side_a, side_b
    elif best.direction == "lower":
        higher_codes, lower_codes = side_b, side_a
    else:
        name_a = ",".join(sorted(family_names[code] for code in side_a))
        name_b = ",".join(sorted(family_names[code] for code in side_b))
        if name_a <= name_b:
            lower_codes, higher_codes = side_a, side_b
        else:
            lower_codes, higher_codes = side_b, side_a
    order = {code: index for index, code in enumerate(range(len(family_names)))}
    lower_codes = sorted(lower_codes, key=lambda code: order[code])
    higher_codes = sorted(higher_codes, key=lambda code: order[code])
    column = sub_values[:, best.column]
    lower_mask = np.isin(sub_groups, np.asarray(lower_codes, dtype=np.int32))
    median_lower = float(np.median(column[lower_mask]))
    median_higher = float(np.median(column[~lower_mask]))
    n_lower = int(sum(int(counts[code]) for code in lower_codes))
    n_higher = int(sum(int(counts[code]) for code in higher_codes))
    expand_lower = best.oriented >= _SPLIT_MIN and len(lower_codes) >= 2
    expand_higher = best.oriented >= _SPLIT_MIN and len(higher_codes) >= 2
    if expand_lower and expand_higher:
        expanded = "both"
    elif expand_lower:
        expanded = "lower"
    elif expand_higher:
        expanded = "higher"
    else:
        expanded = "neither"
    lower_names = tuple(sorted(family_names[code] for code in lower_codes))
    higher_names = tuple(sorted(family_names[code] for code in higher_codes))
    parent = "+".join(sorted(family_names[code] for code in codes))
    split = _Split(
        view=view,
        depth=depth,
        parent=parent,
        feature=best.name,
        oriented_auc=best.oriented,
        threshold=0.5 * (median_lower + median_higher),
        lower_families=lower_names,
        higher_families=higher_names,
        n_lower=n_lower,
        n_higher=n_higher,
        median_lower=median_lower,
        median_higher=median_higher,
        unstable=n_lower < _UNSTABLE_N or n_higher < _UNSTABLE_N,
        expanded=expanded,
    )
    return split, lower_codes, higher_codes


def _subset_tree(
    view: str,
    names: list[str],
    matrix: np.ndarray,
    groups: np.ndarray,
    family_names: list[str],
    counts: np.ndarray,
    columns: list[int],
    epsilon2: np.ndarray,
) -> list[_Split]:
    splits: list[_Split] = []

    def walk(codes: list[int], depth: int) -> None:
        if len(codes) < 2:
            return
        if depth > len(family_names):
            _die("ERROR: subset tree did not shrink")
        split, lower_codes, higher_codes = _best_split(
            codes,
            view,
            depth,
            names,
            matrix,
            groups,
            family_names,
            counts,
            columns,
            epsilon2,
        )
        splits.append(split)
        if split.expanded in ("lower", "both"):
            walk(lower_codes, depth + 1)
        if split.expanded in ("higher", "both"):
            walk(higher_codes, depth + 1)

    walk(list(range(len(family_names))), 1)
    return splits


def _transform(values: np.ndarray, names: list[str]) -> np.ndarray:
    transformed = np.empty_like(values)
    for column, name in enumerate(names):
        column_values = values[:, column]
        if _transform_of(name) == "identity":
            transformed[:, column] = column_values
            continue
        if np.any(column_values < 0.0):
            _die("ERROR: %s is negative and cannot be log1p" % name)
        transformed[:, column] = np.log1p(column_values)
    return transformed


def _folds(groups: np.ndarray, n_groups: int, n_folds: int) -> np.ndarray:
    folds = np.empty(groups.size, dtype=np.int32)
    rng = np.random.RandomState(_SEED)
    for group in range(n_groups):
        index = np.flatnonzero(groups == group)
        rng.shuffle(index)
        for position, row in enumerate(index):
            folds[row] = position % n_folds
    return folds


def _zscore_apply(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = train.mean(axis=0)
    scale = train.std(axis=0, ddof=0)
    zero = scale == 0.0
    scale = scale.copy()
    scale[zero] = 1.0
    train_z = (train - mean) / scale
    test_z = (test - mean) / scale
    if np.any(zero):
        train_z[:, zero] = 0.0
        test_z[:, zero] = 0.0
    return train_z, test_z


def _fit_predict(
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
    n_classes: int,
) -> np.ndarray:
    train_z, test_z = _zscore_apply(train_x, test_x)
    n_rows, n_features = train_z.shape
    design = np.column_stack([train_z, np.ones(n_rows)])
    penalty = np.eye(n_features + 1, dtype=np.float64) * _RIDGE
    penalty[-1, -1] = 0.0
    gram = design.T @ design + penalty
    test_design = np.column_stack([test_z, np.ones(test_x.shape[0])])
    scores = np.empty((test_x.shape[0], n_classes), dtype=np.float64)
    for class_index in range(n_classes):
        target = (train_y == class_index).astype(np.float64)
        try:
            beta = np.linalg.solve(gram, design.T @ target)
        except np.linalg.LinAlgError:
            _die("ERROR: ridge system is singular")
        if not np.all(np.isfinite(beta)):
            _die("ERROR: ridge coefficients are not finite")
        scores[:, class_index] = test_design @ beta
    if not np.all(np.isfinite(scores)):
        _die("ERROR: ridge scores are not finite")
    return scores.argmax(axis=1).astype(np.int32)


def _cross_validate(
    x: np.ndarray,
    y: np.ndarray,
    log_length: np.ndarray,
    folds: np.ndarray,
    n_classes: int,
    residualize: bool,
) -> np.ndarray:
    confusion = np.zeros((n_classes, n_classes), dtype=np.int64)
    for fold in range(int(folds.max()) + 1):
        test = folds == fold
        train = ~test
        if not np.any(test) or not np.any(train):
            _die("ERROR: a CV fold is empty")
        x_train = x[train]
        x_test = x[test]
        if residualize:
            x_train, x_test = _residualize_fold(
                x_train, x_test, log_length[train], log_length[test]
            )
        predicted = _fit_predict(x_train, y[train], x_test, n_classes)
        np.add.at(confusion, (y[test], predicted), 1)
    return confusion


def _balanced_accuracy(confusion: np.ndarray) -> tuple[float, list[float]]:
    recalls: list[float] = []
    for row in confusion:
        support = int(row.sum())
        if support == 0:
            _die("ERROR: a family is missing from the confusion matrix")
        recalls.append(float(row[len(recalls)]) / float(support))
    return float(np.mean(np.asarray(recalls, dtype=np.float64))), recalls


def _pca(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    _train, z = _zscore_apply(x, x)
    _u, singular, _vt = np.linalg.svd(z, full_matrices=False)
    total = float(np.dot(singular, singular))
    if total <= 0.0:
        _die("ERROR: PCA matrix has no variance")
    explained = (singular * singular) / total
    scores = _u * singular
    return scores, explained


def _check_results(
    splits: list[_Split],
    pairs: list[_Pair],
    ridge: list[_Ridge],
    family_names: list[str],
    counts: np.ndarray,
) -> None:
    expected_pairs = len(family_names) * (len(family_names) - 1) // 2
    _expect(len(pairs) == expected_pairs, "pairwise table covers every family pair")
    for pair in pairs:
        _expect(
            _CLEAR_AUC - 1.0 <= pair.oriented_auc <= 1.0 + 1e-9,
            "pairwise AUC is oriented",
        )
        _expect(pair.oriented_auc + 1e-9 >= 0.5, "pairwise oriented AUC is at least 0.5")
        _expect(pair.oriented_auc_adj + 1e-9 >= 0.5, "adjusted pairwise AUC is at least 0.5")
    for split in splits:
        parent = set(split.parent.split("+"))
        lower = set(split.lower_families)
        higher = set(split.higher_families)
        _expect(lower | higher == parent, "split covers its parent")
        _expect(lower.isdisjoint(higher), "split sides are disjoint")
        _expect(0.5 - 1e-9 <= split.oriented_auc <= 1.0 + 1e-9, "split AUC range")
        _expect(
            split.n_lower + split.n_higher == int(sum(counts[family_names.index(name)] for name in parent)),
            "split counts add up",
        )
    for result in ridge:
        for index, name in enumerate(family_names):
            _expect(
                int(result.confusion[index].sum()) == int(counts[index]),
                "confusion support matches %s" % name,
            )


def _shortlist(splits: list[_Split], canon: dict[str, str]) -> list[tuple[str, float, str, str, str]]:
    best: dict[str, list] = {}
    for split in splits:
        if split.oriented_auc < _SPLIT_MIN:
            continue
        name = canon.get(split.feature, split.feature)
        current = best.get(name)
        views = split.view
        if current is None or split.oriented_auc > current[0]:
            view_set = set() if current is None else current[4]
            view_set.add(views)
            best[name] = [
                split.oriented_auc,
                ",".join(split.higher_families),
                ",".join(split.lower_families),
                split.view,
                view_set,
            ]
        else:
            current[4].add(split.view)
    rows = []
    for name, item in best.items():
        view_text = ",".join(sorted(item[4]))
        rows.append((name, float(item[0]), view_text, str(item[1]), str(item[2])))
    rows.sort(key=lambda row: (-row[1], row[0]))
    return rows


def _write_table(path: Path, header: str, rows: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(header)
        if not header.endswith("\n"):
            handle.write("\n")
        for row in rows:
            handle.write(row)
            handle.write("\n")


def _summary(
    n_vle: int,
    n_decoy: int,
    family_names: list[str],
    counts: np.ndarray,
    near_constant: list[tuple[str, float, float]],
    duplicates: list[tuple[str, ...]],
    removed: list[str],
    scores: list[_FeatureScore],
    one_vs_rest: list[_Ovr],
    pairs: list[_Pair],
    splits: list[_Split],
    shortlist: list[tuple[str, float, str, str, str]],
    ridge: list[_Ridge],
    pca_raw: float,
    pca_raw_2: float,
    pca_adj: float,
    pca_adj_2: float,
) -> str:
    lines = [
        "VLE family features. DECOY rows are discarded.",
        "The family is the second field of a VLE_family_... label.",
        "Oriented AUC is max(AUC, 1-AUC). Chance is 0.5.",
        "higher means that family, or that side, tends to have the larger value.",
        "Epsilon-squared is the tie-corrected Kruskal-Wallis H / (n - 1).",
        "Length adjustment residualizes each raw feature except %s on log(%s)."
        % (_LENGTH, _LENGTH),
        "A pair clears %.2f when its best single-feature oriented AUC reaches that value."
        % _CLEAR_AUC,
        "A subset split is expanded when its oriented AUC is at least %.2f and the side has two or more families."
        % _SPLIT_MIN,
        "The threshold is the midpoint of the two side medians.",
        "Sides are named by the Mann-Whitney direction. A tie is ordered alphabetically.",
        "A family or a side with fewer than %d circles is marked unstable." % _UNSTABLE_N,
        "In a pair, family_a is the larger family.",
        "The linear model and the PCA use log1p on counts and lengths, and leave percentages and the named fractions unchanged.",
        "Duplicate copies and near-constant columns are removed from the linear model and the PCA.",
        "Near-constant columns are omitted from length-adjusted scores.",
        "The ridge penalty is %.3g on z-scored predictors. The intercept is not penalized. The fit is unweighted."
        % _RIDGE,
        "Cross-validation is %d-fold, stratified by family, shuffle seed %d."
        % (_FOLDS, _SEED),
        "",
        "classes",
        "VLE\t%d" % n_vle,
        "DECOY_discarded\t%d" % n_decoy,
        "families\t%d" % len(family_names),
    ]
    for index, name in enumerate(family_names):
        lines.append(
            "%s\t%d\tfraction=%s\tunstable=%s"
            % (
                name,
                int(counts[index]),
                _fmt(float(counts[index]) / float(n_vle)),
                _yes(int(counts[index]) < _UNSTABLE_N),
            )
        )
    lines.append("")
    lines.append(
        "near-constant features (at least %.0f%% of VLE rows share one value)"
        % (100.0 * _NEAR_CONSTANT_FRACTION)
    )
    if not near_constant:
        lines.append("none")
    else:
        for name, mode, fraction in near_constant:
            lines.append("%s\tvalue=%s\tfraction=%.4f" % (name, _fmt(mode), fraction))
    lines.append("")
    lines.append("exact duplicate columns")
    if not duplicates:
        lines.append("none")
    else:
        for group in duplicates:
            lines.append("\t".join(group))
    lines.append("")
    lines.append("removed from the linear model and the PCA")
    if not removed:
        lines.append("none")
    else:
        for name in removed:
            lines.append(name)
    lines.append("")
    lines.append("top %d by Kruskal-Wallis epsilon-squared" % _TOP_N)
    for rank, score in enumerate(scores[:_TOP_N], start=1):
        lines.append(
            "%d\t%s\tepsilon2=%.4f\tmax_ovr_auc=%.4f\textreme=%s\tdirection=%s"
            % (
                rank,
                score.name,
                score.epsilon2,
                score.max_ovr,
                score.extreme_family,
                score.extreme_direction,
            )
        )
    adjusted = [score for score in scores if score.epsilon2_adj is not None]
    adjusted.sort(key=lambda score: (-(score.epsilon2_adj or 0.0), score.name))
    lines.append("")
    lines.append(
        "top %d by epsilon-squared after residualizing on log(%s)" % (_TOP_N, _LENGTH)
    )
    lines.append("%s is the covariate and is omitted." % _LENGTH)
    for rank, score in enumerate(adjusted[:_TOP_N], start=1):
        lines.append(
            "%d\t%s\tepsilon2=%.4f\tmax_ovr_auc=%.4f\textreme=%s\tdirection=%s"
            % (
                rank,
                score.name,
                score.epsilon2_adj,
                score.max_ovr_adj,
                score.extreme_family_adj,
                score.extreme_direction_adj,
            )
        )
    lines.append("")
    lines.append("one versus rest, best feature")
    for row in one_vs_rest:
        lines.append(
            "%s\tn=%d\traw=%s\tauc=%.4f\t%s\tadjusted=%s\tauc=%.4f\t%s\tunstable=%s"
            % (
                row.family,
                row.n,
                row.feature,
                row.oriented_auc,
                row.direction,
                row.feature_adj,
                row.oriented_auc_adj,
                row.direction_adj,
                _yes(row.unstable),
            )
        )
    lines.append("")
    lines.append("subset tree, raw")
    lines.extend(_tree_lines(splits, "raw"))
    lines.append("")
    lines.append("subset tree, length-adjusted")
    lines.extend(_tree_lines(splits, "length_adjusted"))
    lines.append("")
    clear = [pair for pair in pairs if pair.oriented_auc >= _CLEAR_AUC]
    mixed = [pair for pair in pairs if pair.oriented_auc < _CLEAR_AUC]
    lines.append("pairs with best raw oriented AUC >= %.2f" % _CLEAR_AUC)
    lines.extend(_pair_lines(clear))
    lines.append("")
    lines.append("pairs that stay mixed on the best raw feature")
    lines.extend(_pair_lines(mixed))
    adjusted_clear = [pair for pair in pairs if pair.oriented_auc_adj >= _CLEAR_AUC]
    raw_only = [
        pair
        for pair in pairs
        if pair.oriented_auc >= _CLEAR_AUC and pair.oriented_auc_adj < _CLEAR_AUC
    ]
    adjusted_only = [
        pair
        for pair in pairs
        if pair.oriented_auc < _CLEAR_AUC and pair.oriented_auc_adj >= _CLEAR_AUC
    ]
    lines.append("")
    lines.append(
        "length: %d of %d pairs clear %.2f on a raw feature, %d of %d still do after residualizing on log(%s), %d clear only before adjustment, %d clear only after adjustment."
        % (
            len(clear),
            len(pairs),
            _CLEAR_AUC,
            len(adjusted_clear),
            len(pairs),
            _LENGTH,
            len(raw_only),
            len(adjusted_only),
        )
    )
    n_length_wins = sum(1 for pair in pairs if pair.feature == _LENGTH)
    lines.append("%s is the best raw feature for %d pairs." % (_LENGTH, n_length_wins))
    lines.append("")
    lines.append("shortlist (features that win a split of oriented AUC >= %.2f)" % _SPLIT_MIN)
    if not shortlist:
        lines.append("none")
    else:
        for name, auc, views, higher, lower in shortlist:
            lines.append(
                "%s\tbest_auc=%.4f\tviews=%s\thigher=%s\tlower=%s"
                % (name, auc, views, higher, lower)
            )
    lines.append("")
    lines.append("ridge one-versus-rest, balanced accuracy")
    by_name = {result.name: result for result in ridge}
    for name in ("nr_bases", "all_features", "length_residualized"):
        result = by_name[name]
        lines.append("%s\t%.4f" % (name, result.balanced_accuracy))
    for name in ("nr_bases", "all_features", "length_residualized"):
        result = by_name[name]
        lines.append("")
        lines.append("per-family recall, %s" % name)
        for index, family in enumerate(family_names):
            lines.append("%s\t%.4f" % (family, result.recalls[index]))
    lines.append("")
    lines.append("pca variance explained, percent")
    lines.append("raw_transformed\tpc1=%.1f\tpc2=%.1f" % (100.0 * pca_raw, 100.0 * pca_raw_2))
    lines.append(
        "length_adjusted\tpc1=%.1f\tpc2=%.1f" % (100.0 * pca_adj, 100.0 * pca_adj_2)
    )
    lines.append("PCA points are capped at %d circles per family. Centroids use every circle." % _PCA_PER_FAMILY)
    lines.append("")
    lines.extend(
        _next_steps(
            pairs,
            splits,
            by_name,
            shortlist,
            one_vs_rest,
            n_length_wins,
            len(raw_only),
        )
    )
    lines.append("")
    lines.append("Full tables: family_counts.tsv, feature_family_scores.tsv, pairwise_best_auc.tsv, one_vs_rest.tsv, subset_splits.tsv, family_medians.tsv, confusion_ridge.tsv")
    lines.append(
        "Figures: feature_effect.png, pairwise_best_auc.png, family_median_heatmap.png, pca_families.png, pca_length_adjusted.png, box_top_features.png, box_residual_features.png, box_uncorrelated_features.png, confusion_ridge.png"
    )
    lines.append("")
    return "\n".join(lines)


def _tree_lines(splits: list[_Split], view: str) -> list[str]:
    rows = [split for split in splits if split.view == view]
    if not rows:
        return ["none"]
    lines = []
    for split in rows:
        indent = "  " * (split.depth - 1)
        lines.append(
            "%sdepth=%d\tauc=%.4f\tfeature=%s\tthreshold=%s\thigher=%s\tlower=%s\tn_higher=%d\tn_lower=%d\tmedian_higher=%s\tmedian_lower=%s\texpanded=%s\tunstable=%s"
            % (
                indent,
                split.depth,
                split.oriented_auc,
                split.feature,
                _fmt(split.threshold),
                ",".join(split.higher_families),
                ",".join(split.lower_families),
                split.n_higher,
                split.n_lower,
                _fmt(split.median_higher),
                _fmt(split.median_lower),
                split.expanded,
                _yes(split.unstable),
            )
        )
    return lines


def _pair_lines(pairs: list[_Pair]) -> list[str]:
    if not pairs:
        return ["none"]
    lines = []
    for pair in pairs:
        lines.append(
            "%s\t%s\tn=%d+%d\traw=%s\tauc=%.4f\tfamily_a_values=%s\tnr_bases_auc=%.4f\tadjusted=%s\tauc=%.4f\tfamily_a_values_adjusted=%s\tunstable=%s"
            % (
                pair.family_a,
                pair.family_b,
                pair.n_a,
                pair.n_b,
                pair.feature,
                pair.oriented_auc,
                pair.family_a_direction,
                pair.nr_bases_oriented_auc,
                pair.feature_adj,
                pair.oriented_auc_adj,
                pair.family_a_direction_adj,
                _yes(pair.unstable),
            )
        )
    return lines


def _next_steps(
    pairs: list[_Pair],
    splits: list[_Split],
    ridge: dict[str, _Ridge],
    shortlist: list[tuple[str, float, str, str, str]],
    one_vs_rest: list[_Ovr],
    n_length_wins: int,
    n_raw_only: int,
) -> list[str]:
    lines = ["next classifier"]
    bal_length = ridge["nr_bases"].balanced_accuracy
    bal_all = ridge["all_features"].balanced_accuracy
    bal_adj = ridge["length_residualized"].balanced_accuracy
    lines.append(
        "Balanced accuracy is %.3f from %s alone, %.3f from all kept features, and %.3f after removing log length."
        % (bal_length, _LENGTH, bal_all, bal_adj)
    )
    if bal_all >= bal_length + 0.05:
        lines.append("The other features improve on length in the linear model.")
    else:
        lines.append("The full linear model is close to length alone, so size is the main linear signal.")
    adjusted_winners = [pair for pair in pairs if pair.oriented_auc_adj >= _CLEAR_AUC]
    count_winners = [
        pair for pair in adjusted_winners if _transform_of(pair.feature_adj) == "log1p"
    ]
    fraction_winners = [
        pair for pair in adjusted_winners if _transform_of(pair.feature_adj) != "log1p"
    ]
    lines.append(
        "After residualizing on log(%s), %d of %d pairs still reach %.2f."
        % (_LENGTH, len(adjusted_winners), len(pairs), _CLEAR_AUC)
    )
    lines.append(
        "%d of those winners are residual counts. A count residual can still follow circle size."
        % len(count_winners)
    )
    if not fraction_winners:
        lines.append(
            "No percentage or fraction reaches %.2f on a length residual." % _CLEAR_AUC
        )
    elif len(fraction_winners) == 1:
        pair = fraction_winners[0]
        lines.append(
            "One winner is a percentage or fraction: %s versus %s on %s (%.3f)."
            % (pair.family_a, pair.family_b, pair.feature_adj, pair.oriented_auc_adj)
        )
    else:
        described = ", ".join(
            "%s versus %s on %s (%.3f)"
            % (pair.family_a, pair.family_b, pair.feature_adj, pair.oriented_auc_adj)
            for pair in fraction_winners
        )
        lines.append(
            "%d winners are percentages or fractions: %s."
            % (len(fraction_winners), described)
        )
    ovr_names = [row.feature for row in one_vs_rest]
    if ovr_names and len(set(ovr_names)) == 1:
        lines.append(
            "%s is the best raw one-versus-rest feature for every family. Try it in a family-versus-rest classifier even when the subset tree splits on another feature first."
            % ovr_names[0]
        )
    if n_raw_only > 0:
        lines.append(
            "%d pairs reach %.2f only before length adjustment."
            % (n_raw_only, _CLEAR_AUC)
        )
    if n_length_wins > 0:
        lines.append(
            "Where %s wins the raw comparison, put an explicit length term in the classifier or use the length-adjusted tree."
            % _LENGTH
        )
    raw_root = [split for split in splits if split.view == "raw" and split.depth == 1]
    adj_root = [split for split in splits if split.view == "length_adjusted" and split.depth == 1]
    if raw_root:
        split = raw_root[0]
        lines.append(
            "The best raw split is oriented AUC %.3f on %s: higher %s versus lower %s."
            % (
                split.oriented_auc,
                split.feature,
                ",".join(split.higher_families),
                ",".join(split.lower_families),
            )
        )
    if adj_root:
        split = adj_root[0]
        lines.append(
            "The best length-adjusted split is oriented AUC %.3f on %s: higher %s versus lower %s."
            % (
                split.oriented_auc,
                split.feature,
                ",".join(split.higher_families),
                ",".join(split.lower_families),
            )
        )
    if shortlist:
        lines.append(
            "Start with the shortlist rather than all columns: %s."
            % ", ".join(row[0] for row in shortlist)
        )
    else:
        lines.append(
            "No split reached %.2f, so there is no single-feature partition to start from."
            % _SPLIT_MIN
        )
    lines.append(
        "The ridge fit is unweighted, so the largest family dominates the loss. Read balanced accuracy and the subset tree before raw accuracy."
    )
    lines.append(
        "zetavirus has fewer than %d circles. Keep it in the tables, and do not treat a zetavirus call as established."
        % _UNSTABLE_N
    )
    lines.append(
        "Use the confusion matrix to see which families still mix. A classifier of two or more subsets is the goal if the off-block families are the mixed ones."
    )
    return lines


def _save_fig(fig: plt.Figure, path: Path) -> None:
    temporary = path.with_suffix(".tmp.png")
    fig.savefig(temporary, dpi=150, bbox_inches="tight")
    plt.close(fig)
    temporary.replace(path)


def _write_feature_effect(path: Path, scores: list[_FeatureScore]) -> None:
    ordered = sorted(
        scores,
        key=lambda score: (max(score.epsilon2, score.epsilon2_adj or 0.0), score.name),
    )
    names = [score.name for score in ordered]
    raw = np.asarray([score.epsilon2 for score in ordered], dtype=np.float64)
    adjusted = np.asarray(
        [np.nan if score.epsilon2_adj is None else score.epsilon2_adj for score in ordered],
        dtype=np.float64,
    )
    positions = np.arange(len(names))
    fig_height = max(6.0, 0.36 * len(names) + 1.5)
    fig, ax = plt.subplots(figsize=(10.5, fig_height))
    ax.barh(positions + 0.16, raw, height=0.32, color=_COLOR_RAW, label="raw", zorder=2)
    present = np.isfinite(adjusted)
    ax.barh(
        positions[present] - 0.16,
        adjusted[present],
        height=0.32,
        color=_COLOR_ADJ,
        label="length-adjusted",
        zorder=2,
    )
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("Epsilon-squared")
    ax.set_title("Kruskal-Wallis epsilon-squared across VLE families")
    ax.grid(axis="x", color="0.85", zorder=0)
    ax.tick_params(axis="y", length=0)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    _save_fig(fig, path)


def _write_pairwise_heatmap(
    path: Path, pairs: list[_Pair], family_names: list[str]
) -> None:
    n_families = len(family_names)
    index = {name: position for position, name in enumerate(family_names)}
    raw = np.full((n_families, n_families), np.nan, dtype=np.float64)
    adjusted = np.full((n_families, n_families), np.nan, dtype=np.float64)
    for pair in pairs:
        i = index[pair.family_a]
        j = index[pair.family_b]
        raw[i, j] = pair.oriented_auc
        raw[j, i] = pair.oriented_auc
        adjusted[i, j] = pair.oriented_auc_adj
        adjusted[j, i] = pair.oriented_auc_adj
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.8))
    cmap = plt.cm.viridis.copy()
    cmap.set_bad("0.9")
    image = None
    for ax, matrix, title in (
        (axes[0], raw, "Best raw feature"),
        (axes[1], adjusted, "Best length-adjusted feature"),
    ):
        image = ax.imshow(
            np.ma.masked_invalid(matrix),
            vmin=0.5,
            vmax=1.0,
            cmap=cmap,
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_xticks(np.arange(n_families))
        ax.set_yticks(np.arange(n_families))
        ax.set_xticklabels(family_names, rotation=30, ha="right", fontsize=8)
        ax.set_yticklabels(family_names, fontsize=8)
        ax.set_title(title)
        for i in range(n_families):
            for j in range(n_families):
                if i == j or not math.isfinite(matrix[i, j]):
                    continue
                value = float(matrix[i, j])
                color = "white" if value < 0.75 else "black"
                ax.text(j, i, "%.2f" % value, ha="center", va="center", fontsize=7, color=color)
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes, fraction=0.046, pad=0.04)
        colorbar.set_label("Oriented AUC")
    fig.suptitle("Best single-feature separation of each family pair", fontsize=13)
    fig.subplots_adjust(top=0.86, wspace=0.55)
    _save_fig(fig, path)


def _family_z(medians: np.ndarray) -> np.ndarray:
    z = np.empty_like(medians)
    for column in range(medians.shape[1]):
        values = medians[:, column]
        scale = float(values.std(ddof=0))
        if scale == 0.0:
            z[:, column] = 0.0
        else:
            z[:, column] = (values - float(values.mean())) / scale
    return z


def _cluster_order(z: np.ndarray, columns: list[int]) -> list[int]:
    if not columns:
        return list(range(z.shape[0]))
    profile = z[:, columns]
    difference = profile[:, None, :] - profile[None, :, :]
    distance = np.sqrt(np.sum(difference * difference, axis=2))
    return _leaf_order(distance)


def _write_median_heatmap(
    path: Path,
    medians: np.ndarray,
    family_names: list[str],
    feature_names: list[str],
    scores: list[_FeatureScore],
    cluster_columns: list[int],
) -> None:
    z = _family_z(medians)
    order = _cluster_order(z, cluster_columns)
    feature_order = sorted(
        range(len(feature_names)),
        key=lambda index: (
            -next(score.epsilon2 for score in scores if score.name == feature_names[index]),
            feature_names[index],
        ),
    )
    matrix = z[np.ix_(order, feature_order)].T
    limit = float(np.max(np.abs(matrix))) if matrix.size else 1.0
    if limit == 0.0:
        limit = 1.0
    fig_height = max(6.0, 0.28 * len(feature_names) + 1.8)
    fig, ax = plt.subplots(figsize=(8.2, fig_height))
    image = ax.imshow(
        matrix,
        vmin=-limit,
        vmax=limit,
        cmap="coolwarm",
        interpolation="nearest",
        aspect="auto",
    )
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([family_names[index] for index in order], rotation=30, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(feature_order)))
    ax.set_yticklabels([feature_names[index] for index in feature_order], fontsize=7)
    ax.set_title("Family medians, z-scored across families")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Z-scored median")
    fig.tight_layout()
    _save_fig(fig, path)


def _write_pca(
    path: Path,
    scores: np.ndarray,
    explained: np.ndarray,
    groups: np.ndarray,
    family_names: list[str],
    title: str,
) -> None:
    if scores.shape[1] < 2:
        _die("ERROR: PCA has fewer than two components")
    rng = np.random.RandomState(_SEED)
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    for group, name in enumerate(family_names):
        index = np.flatnonzero(groups == group)
        if index.size > _PCA_PER_FAMILY:
            index = rng.choice(index, size=_PCA_PER_FAMILY, replace=False)
        color = _COLORS[group % len(_COLORS)]
        ax.scatter(
            scores[index, 0],
            scores[index, 1],
            s=12,
            alpha=0.45,
            color=color,
            linewidths=0,
            label=name,
            zorder=2,
        )
        centroid = scores[groups == group, :2].mean(axis=0)
        ax.scatter(
            [centroid[0]],
            [centroid[1]],
            s=70,
            color=color,
            edgecolors="black",
            linewidths=0.6,
            zorder=3,
        )
        ax.annotate(
            name,
            centroid,
            textcoords="offset points",
            xytext=(5, 4),
            fontsize=8,
            zorder=4,
        )
    ax.set_xlabel("PC1 (%.1f%%)" % (100.0 * float(explained[0])))
    ax.set_ylabel("PC2 (%.1f%%)" % (100.0 * float(explained[1])))
    ax.set_title(title)
    ax.axhline(0.0, color="0.85", linewidth=0.8, zorder=0)
    ax.axvline(0.0, color="0.85", linewidth=0.8, zorder=0)
    ax.legend(frameon=True, fontsize=8, loc="best")
    fig.tight_layout()
    _save_fig(fig, path)


def _length_correlation(values: np.ndarray, length_index: int) -> np.ndarray:
    """Spearman rho of every column with nr_bases."""
    length_ranks = _average_ranks(values[:, length_index])
    centered_length = length_ranks - length_ranks.mean()
    length_norm = float(np.sqrt(np.dot(centered_length, centered_length)))
    rho = np.empty(values.shape[1], dtype=np.float64)
    for index in range(values.shape[1]):
        ranks = _average_ranks(values[:, index])
        centered = ranks - ranks.mean()
        column_norm = float(np.sqrt(np.dot(centered, centered)))
        if length_norm == 0.0 or column_norm == 0.0:
            rho[index] = 0.0
        else:
            rho[index] = float(np.dot(centered_length, centered) / (length_norm * column_norm))
    return rho


def _uncorrelated_features(
    scores: list[_FeatureScore],
    names: list[str],
    rho: np.ndarray,
    near_names: set[str],
    canon: dict[str, str],
) -> list[tuple[_FeatureScore, float]]:
    by_name = {score.name: score for score in scores}
    chosen: list[tuple[float, _FeatureScore, float]] = []
    for index, name in enumerate(names):
        if name == _LENGTH or name in near_names:
            continue
        if canon.get(name, name) != name:
            continue
        correlation = float(rho[index])
        if abs(correlation) >= _LENGTH_RHO_MAX:
            continue
        chosen.append((abs(correlation), by_name[name], correlation))
    chosen.sort(key=lambda item: (item[0], item[1].name))
    return [(score, correlation) for _abs_rho, score, correlation in chosen]


def _write_boxes(
    path: Path,
    values: np.ndarray,
    groups: np.ndarray,
    family_names: list[str],
    names: list[str],
    chosen: list[_FeatureScore],
    title: str,
    residual: bool,
    notes: list[str] | None = None,
) -> None:
    if not chosen:
        _die("ERROR: no features to boxplot")
    ncols = 3 if len(chosen) > 3 else len(chosen)
    nrows = int(math.ceil(len(chosen) / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.1 * ncols, 3.3 * nrows), squeeze=False)
    flat = axes.ravel()
    name_to_index = {name: index for index, name in enumerate(names)}
    positions = np.arange(1, len(family_names) + 1)
    for panel, score in enumerate(chosen):
        ax = flat[panel]
        epsilon = score.epsilon2_adj if residual else score.epsilon2
        if epsilon is None:
            _die("ERROR: %s has no epsilon-squared for this box plot" % score.name)
        column = values[:, name_to_index[score.name]]
        data = [column[groups == group] for group in range(len(family_names))]
        drawn = ax.boxplot(data, positions=positions, showfliers=False, patch_artist=True)
        for box, group in zip(drawn["boxes"], range(len(family_names))):
            box.set_facecolor(_COLORS[group % len(_COLORS)])
            box.set_alpha(0.8)
        if residual:
            ax.axhline(0.0, color="0.35", linewidth=0.8, zorder=0)
        ax.set_xticks(positions)
        ax.set_xticklabels(family_names, rotation=30, ha="right", fontsize=7)
        if notes is not None:
            ax.set_title(notes[panel], fontsize=8)
        else:
            ax.set_title("%s\nepsilon2 %.3f" % (score.name, epsilon), fontsize=9)
        ax.tick_params(labelsize=7)
        if panel % ncols == 0:
            ax.set_ylabel("Residual" if residual else "Value", fontsize=8)
    for panel in range(len(chosen), flat.size):
        flat[panel].set_visible(False)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    _save_fig(fig, path)


def _write_confusion(
    path: Path,
    confusion: np.ndarray,
    family_names: list[str],
    balanced_accuracy: float,
) -> None:
    support = confusion.sum(axis=1, keepdims=True).astype(np.float64)
    normalized = np.zeros_like(confusion, dtype=np.float64)
    np.divide(confusion, support, out=normalized, where=support > 0)
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    image = ax.imshow(normalized, vmin=0.0, vmax=1.0, cmap="Blues", interpolation="nearest")
    ax.set_xticks(np.arange(len(family_names)))
    ax.set_yticks(np.arange(len(family_names)))
    ax.set_xticklabels(family_names, rotation=30, ha="right", fontsize=8)
    ax.set_yticklabels(family_names, fontsize=8)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title("All features, row-normalized\nbalanced accuracy %.3f" % balanced_accuracy)
    for i in range(len(family_names)):
        for j in range(len(family_names)):
            value = float(normalized[i, j])
            color = "white" if value >= 0.55 else "black"
            ax.text(j, i, "%.2f" % value, ha="center", va="center", fontsize=8, color=color)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Recall fraction")
    fig.tight_layout()
    _save_fig(fig, path)


def _box_features(
    scores: list[_FeatureScore],
    near_names: set[str],
    canon: dict[str, str],
    residual: bool,
) -> list[_FeatureScore]:
    if residual:
        ordered = [score for score in scores if score.epsilon2_adj is not None]
        ordered.sort(key=lambda score: (-(score.epsilon2_adj or 0.0), score.name))
    else:
        ordered = scores
    chosen: list[_FeatureScore] = []
    seen: set[str] = set()
    for score in ordered:
        if score.name in near_names:
            continue
        if residual and score.name == _LENGTH:
            continue
        name = canon.get(score.name, score.name)
        if name in seen:
            continue
        seen.add(name)
        chosen.append(score)
        if len(chosen) == _BOX_N:
            break
    return chosen


def _self_check() -> None:
    _expect(_family_of("VLE_ambivirus_rest", 1) == "ambivirus", "family is the second field")
    _expect(_family_of("DECOY_ambivirus_rest", 1) is None, "DECOY rows are discarded")
    _expect(_transform_of("pct_paired_bases") == "identity", "percentage stays raw")
    _expect(_transform_of("pRod") == "identity", "pRod stays raw")
    _expect(_transform_of("nBranch") == "log1p", "branch count is log1p")
    _expect(_transform_of(_LENGTH) == "log1p", "length is log1p")
    perfect = _auc_from_ranks(
        _average_ranks(np.array([1.0, 2.0, 3.0, 4.0])),
        np.array([False, False, True, True]),
    )
    _expect(abs(perfect - 1.0) < 1e-12, "perfect AUC")
    tied = _auc_from_ranks(
        _average_ranks(np.array([1.0, 1.0, 1.0, 1.0])),
        np.array([True, False, True, False]),
    )
    _expect(abs(tied - 0.5) < 1e-12, "tied AUC")
    constant_h, constant_eps = _kruskal(
        np.array([5.0, 5.0, 5.0, 5.0]),
        np.array([0, 0, 1, 1]),
        2,
    )
    _expect(constant_h == 0.0 and constant_eps == 0.0, "constant Kruskal-Wallis")
    separated_h, separated_eps = _kruskal(
        np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0]),
        np.array([0, 0, 0, 1, 1, 1]),
        2,
    )
    _expect(abs(separated_h - 5.0) < 1e-8, "tie-corrected Kruskal-Wallis")
    _expect(abs(separated_eps - 1.0) < 1e-12, "epsilon-squared of a perfect split")
    _expect(len(_bipartitions([0, 1, 2, 3, 4, 5])) == 31, "six families have 31 splits")
    _expect(len(_bipartitions([0, 1, 2])) == 3, "three families have 3 splits")
    distance = np.array(
        [[0.0, 1.0, 10.0], [1.0, 0.0, 10.0], [10.0, 10.0, 0.0]],
        dtype=np.float64,
    )
    _expect(_leaf_order(distance) == [0, 1, 2], "average linkage keeps neighbors together")
    length = np.array([1.0, 2.0, 4.0, 8.0])
    log_length = np.log(length)
    column = 3.0 + 2.0 * log_length
    residual = _residualize(column.reshape(4, 1), log_length)
    _expect(float(np.max(np.abs(residual))) < 1e-8, "residual of a log-length line is zero")
    toy_x = np.array([[0.0], [0.0], [0.0], [5.0], [5.0], [5.0]])
    toy_y = np.array([0, 0, 0, 1, 1, 1], dtype=np.int32)
    predicted = _fit_predict(toy_x, toy_y, toy_x, 2)
    _expect(np.array_equal(predicted, toy_y), "ridge separates a one-feature gap")
    toy_groups = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1], dtype=np.int32)
    toy_folds = _folds(toy_groups, 2, 5)
    for fold in range(5):
        in_fold = toy_groups[toy_folds == fold]
        _expect(int(np.sum(in_fold == 0)) == 1, "fold keeps one row of class 0")
        _expect(int(np.sum(in_fold == 1)) == 1, "fold keeps one row of class 1")


def main() -> None:
    started = time.time()
    _self_check()
    if len(sys.argv) > 3:
        _die("ERROR: usage: explore_VLE_family_classifier.py [input.tsv] [output_dir]")
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else _DEFAULT_INPUT
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else _DEFAULT_OUTPUT
    feature_names, values, family_list, n_decoy = _load(input_path)
    counts_by_name: dict[str, int] = {}
    for family in family_list:
        counts_by_name[family] = counts_by_name.get(family, 0) + 1
    if len(counts_by_name) < 2:
        _die("ERROR: need at least two VLE families")
    family_names = sorted(counts_by_name, key=lambda name: (-counts_by_name[name], name))
    name_to_code = {name: index for index, name in enumerate(family_names)}
    groups = np.asarray([name_to_code[family] for family in family_list], dtype=np.int32)
    counts = np.bincount(groups, minlength=len(family_names))
    length_index = feature_names.index(_LENGTH)
    length = values[:, length_index]
    if np.any(length <= 0.0):
        _die("ERROR: %s must be positive" % _LENGTH)
    log_length = np.log(length)
    residual = _residualize(values, log_length)
    near_constant = _near_constant(feature_names, values)
    near_names = {item[0] for item in near_constant}
    duplicates = _duplicate_groups(values, feature_names)
    canon = {name: name for name in feature_names}
    for group in duplicates:
        for name in group:
            canon[name] = group[0]
    removed = []
    keep_names = []
    for name in feature_names:
        if name in near_names or canon[name] != name:
            removed.append(name)
            continue
        keep_names.append(name)
    if _LENGTH not in keep_names:
        _die("ERROR: %s was removed from the linear model" % _LENGTH)
    eligible = [index for index, name in enumerate(feature_names) if name not in near_names]
    eligible_adj = [index for index in eligible if index != length_index]
    if not eligible or not eligible_adj:
        _die("ERROR: no varying features to score")
    print(
        "loaded VLE %d DECOY %d families %d"
        % (int(counts.sum()), n_decoy, len(family_names)),
        file=sys.stderr,
        flush=True,
    )
    scores, one_vs_rest = _score_features(
        feature_names,
        values,
        residual,
        groups,
        family_names,
        length_index,
        eligible,
        eligible_adj,
        near_names,
    )
    epsilon2 = np.zeros(len(feature_names), dtype=np.float64)
    epsilon2_adj = np.zeros(len(feature_names), dtype=np.float64)
    for score in scores:
        column = feature_names.index(score.name)
        epsilon2[column] = score.epsilon2
        epsilon2_adj[column] = 0.0 if score.epsilon2_adj is None else score.epsilon2_adj
    pairs = _pairwise(
        feature_names,
        values,
        residual,
        groups,
        family_names,
        counts,
        length_index,
        eligible,
        eligible_adj,
        epsilon2,
        epsilon2_adj,
    )
    raw_splits = _subset_tree(
        "raw",
        feature_names,
        values,
        groups,
        family_names,
        counts,
        eligible,
        epsilon2,
    )
    adjusted_splits = _subset_tree(
        "length_adjusted",
        feature_names,
        residual,
        groups,
        family_names,
        counts,
        eligible_adj,
        epsilon2_adj,
    )
    splits = raw_splits + adjusted_splits
    print(
        "scored %d features, %d pairs, %d splits"
        % (len(scores), len(pairs), len(splits)),
        file=sys.stderr,
        flush=True,
    )
    transformed = _transform(values, feature_names)
    keep_index = [feature_names.index(name) for name in keep_names]
    structure_names = [name for name in keep_names if name != _LENGTH]
    structure_index = [feature_names.index(name) for name in structure_names]
    x_all = transformed[:, keep_index]
    x_length = transformed[:, [length_index]]
    x_structure = transformed[:, structure_index]
    folds = _folds(groups, len(family_names), _FOLDS)
    ridge_specs = (
        ("nr_bases", x_length, False),
        ("all_features", x_all, False),
        ("length_residualized", x_structure, True),
    )
    ridge: list[_Ridge] = []
    for name, matrix, do_residual in ridge_specs:
        confusion = _cross_validate(
            matrix, groups, log_length, folds, len(family_names), do_residual
        )
        balanced, recalls = _balanced_accuracy(confusion)
        ridge.append(_Ridge(name, balanced, recalls, confusion))
        print("ridge %s balanced_accuracy=%.4f" % (name, balanced), file=sys.stderr, flush=True)
    _check_results(splits, pairs, ridge, family_names, counts)
    pca_scores, pca_explained = _pca(x_all)
    structure_residual = _residualize(x_structure, log_length)
    pca_adj_scores, pca_adj_explained = _pca(structure_residual)
    medians = np.vstack(
        [np.median(values[groups == group], axis=0) for group in range(len(family_names))]
    )
    shortlist = _shortlist(splits, canon)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_outputs(
        output_dir,
        family_names,
        counts,
        int(counts.sum()),
        n_decoy,
        feature_names,
        near_constant,
        duplicates,
        removed,
        scores,
        one_vs_rest,
        pairs,
        splits,
        shortlist,
        ridge,
        medians,
        values,
        groups,
        keep_index,
        pca_scores,
        pca_explained,
        pca_adj_scores,
        pca_adj_explained,
        residual,
    )
    elapsed = time.time() - started
    print("finished in %.1fs" % elapsed, file=sys.stderr, flush=True)


def _write_outputs(
    output_dir: Path,
    family_names: list[str],
    counts: np.ndarray,
    n_vle: int,
    n_decoy: int,
    feature_names: list[str],
    near_constant: list[tuple[str, float, float]],
    duplicates: list[tuple[str, ...]],
    removed: list[str],
    scores: list[_FeatureScore],
    one_vs_rest: list[_Ovr],
    pairs: list[_Pair],
    splits: list[_Split],
    shortlist: list[tuple[str, float, str, str, str]],
    ridge: list[_Ridge],
    medians: np.ndarray,
    values: np.ndarray,
    groups: np.ndarray,
    keep_index: list[int],
    pca_scores: np.ndarray,
    pca_explained: np.ndarray,
    pca_adj_scores: np.ndarray,
    pca_adj_explained: np.ndarray,
    residual: np.ndarray,
) -> None:
    n_vle_float = float(n_vle)
    count_rows = [
        "%s\t%d\t%s\t%s"
        % (
            name,
            int(counts[index]),
            _fmt(float(counts[index]) / n_vle_float),
            _yes(int(counts[index]) < _UNSTABLE_N),
        )
        for index, name in enumerate(family_names)
    ]
    _write_table(
        output_dir / "family_counts.tsv",
        "family\tn\tfraction\tunstable\n",
        count_rows,
    )
    score_rows = []
    for score in scores:
        score_rows.append(
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
            % (
                score.name,
                _fmt(score.epsilon2),
                _fmt(score.h),
                _fmt(score.max_ovr),
                score.extreme_family,
                score.extreme_direction,
                _fmt(score.epsilon2_adj),
                _fmt(score.h_adj),
                _fmt(score.max_ovr_adj),
                "NA" if score.extreme_family_adj is None else score.extreme_family_adj,
                "NA" if score.extreme_direction_adj is None else score.extreme_direction_adj,
            )
        )
    _write_table(
        output_dir / "feature_family_scores.tsv",
        "feature\tepsilon2\tH\tmax_ovr_auc\textreme_family\textreme_direction\t"
        "epsilon2_length_adjusted\tH_length_adjusted\tmax_ovr_auc_length_adjusted\t"
        "extreme_family_length_adjusted\textreme_direction_length_adjusted\n",
        score_rows,
    )
    pair_rows = [
        "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
        % (
            pair.family_a,
            pair.family_b,
            pair.n_a,
            pair.n_b,
            _yes(pair.unstable),
            pair.feature,
            _fmt(pair.oriented_auc),
            pair.family_a_direction,
            _fmt(pair.nr_bases_oriented_auc),
            pair.feature_adj,
            _fmt(pair.oriented_auc_adj),
            pair.family_a_direction_adj,
        )
        for pair in pairs
    ]
    _write_table(
        output_dir / "pairwise_best_auc.tsv",
        "family_a\tfamily_b\tn_a\tn_b\tunstable\tfeature\toriented_auc\tfamily_a_values\t"
        "nr_bases_oriented_auc\tfeature_length_adjusted\toriented_auc_length_adjusted\t"
        "family_a_values_length_adjusted\n",
        pair_rows,
    )
    ovr_rows = [
        "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
        % (
            row.family,
            row.n,
            _yes(row.unstable),
            row.feature,
            _fmt(row.oriented_auc),
            row.direction,
            row.feature_adj,
            _fmt(row.oriented_auc_adj),
            row.direction_adj,
        )
        for row in one_vs_rest
    ]
    _write_table(
        output_dir / "one_vs_rest.tsv",
        "family\tn\tunstable\tfeature\toriented_auc\tfamily_values\t"
        "feature_length_adjusted\toriented_auc_length_adjusted\tfamily_values_length_adjusted\n",
        ovr_rows,
    )
    split_rows = [
        "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s"
        % (
            split.view,
            split.depth,
            split.parent,
            split.feature,
            _fmt(split.oriented_auc),
            _fmt(split.threshold),
            ",".join(split.lower_families),
            ",".join(split.higher_families),
            split.n_lower,
            split.n_higher,
            _fmt(split.median_lower),
            _fmt(split.median_higher),
            _yes(split.unstable),
            split.expanded,
        )
        for split in splits
    ]
    _write_table(
        output_dir / "subset_splits.tsv",
        "view\tdepth\tparent\tfeature\toriented_auc\tthreshold\tlower_families\thigher_families\t"
        "n_lower\tn_higher\tmedian_lower\tmedian_higher\tunstable\texpanded\n",
        split_rows,
    )
    median_header = "family\tn\t" + "\t".join(feature_names)
    median_rows = []
    for group, name in enumerate(family_names):
        cells = [name, str(int(counts[group]))]
        cells.extend(_fmt(float(value)) for value in medians[group])
        median_rows.append("\t".join(cells))
    _write_table(output_dir / "family_medians.tsv", median_header + "\n", median_rows)
    confusion_rows = []
    for result in ridge:
        for i, true_name in enumerate(family_names):
            for j, pred_name in enumerate(family_names):
                confusion_rows.append(
                    "%s\t%s\t%s\t%d"
                    % (result.name, true_name, pred_name, int(result.confusion[i, j]))
                )
    _write_table(
        output_dir / "confusion_ridge.tsv",
        "ablation\ttrue_family\tpredicted_family\tcount\n",
        confusion_rows,
    )
    near_set = {item[0] for item in near_constant}
    canon = {name: name for name in feature_names}
    for group in duplicates:
        for name in group:
            canon[name] = group[0]
    length_index = feature_names.index(_LENGTH)
    length_rho = _length_correlation(values, length_index)
    uncorrelated = _uncorrelated_features(scores, feature_names, length_rho, near_set, canon)
    text = _summary(
        n_vle,
        n_decoy,
        family_names,
        counts,
        near_constant,
        duplicates,
        removed,
        scores,
        one_vs_rest,
        pairs,
        splits,
        shortlist,
        ridge,
        float(pca_explained[0]),
        float(pca_explained[1]),
        float(pca_adj_explained[0]),
        float(pca_adj_explained[1]),
    )
    uncorrelated_lines = [
        "features with |Spearman rho| < %.2f versus %s" % (_LENGTH_RHO_MAX, _LENGTH),
        "Duplicate copies and near-constant columns are omitted.",
    ]
    if not uncorrelated:
        uncorrelated_lines.append("none")
    else:
        for score, correlation in uncorrelated:
            uncorrelated_lines.append(
                "%s\trho=%+.3f\tepsilon2=%.4f" % (score.name, correlation, score.epsilon2)
            )
    uncorrelated_lines.append("")
    text = text.replace("Full tables:", "\n".join(uncorrelated_lines) + "Full tables:", 1)
    summary_path = output_dir / "summary.txt"
    summary_path.write_text(text, encoding="utf-8", newline="\n")
    _write_feature_effect(output_dir / "feature_effect.png", scores)
    _write_pairwise_heatmap(output_dir / "pairwise_best_auc.png", pairs, family_names)
    _write_median_heatmap(
        output_dir / "family_median_heatmap.png",
        medians,
        family_names,
        feature_names,
        scores,
        keep_index,
    )
    _write_pca(
        output_dir / "pca_families.png",
        pca_scores,
        pca_explained,
        groups,
        family_names,
        "VLE families, log1p counts, z-scored",
    )
    _write_pca(
        output_dir / "pca_length_adjusted.png",
        pca_adj_scores,
        pca_adj_explained,
        groups,
        family_names,
        "Residualized on log(%s), %s omitted" % (_LENGTH, _LENGTH),
    )
    _write_boxes(
        output_dir / "box_top_features.png",
        values,
        groups,
        family_names,
        feature_names,
        _box_features(scores, near_set, canon, residual=False),
        "Strongest raw features by family",
        residual=False,
    )
    _write_boxes(
        output_dir / "box_residual_features.png",
        residual,
        groups,
        family_names,
        feature_names,
        _box_features(scores, near_set, canon, residual=True),
        "Strongest features after residualizing on log(%s)" % _LENGTH,
        residual=True,
    )
    if uncorrelated:
        _write_boxes(
            output_dir / "box_uncorrelated_features.png",
            values,
            groups,
            family_names,
            feature_names,
            [score for score, _correlation in uncorrelated],
            "Features with |Spearman rho| < %.2f versus %s" % (_LENGTH_RHO_MAX, _LENGTH),
            residual=False,
            notes=[
                "%s\nrho %+.2f   epsilon2 %.3f" % (score.name, correlation, score.epsilon2)
                for score, correlation in uncorrelated
            ],
        )
    all_features = next(result for result in ridge if result.name == "all_features")
    _write_confusion(
        output_dir / "confusion_ridge.png",
        all_features.confusion,
        family_names,
        all_features.balanced_accuracy,
    )
    print(text, end="")
    for name in (
        "summary.txt",
        "family_counts.tsv",
        "feature_family_scores.tsv",
        "pairwise_best_auc.tsv",
        "one_vs_rest.tsv",
        "subset_splits.tsv",
        "family_medians.tsv",
        "confusion_ridge.tsv",
        "feature_effect.png",
        "pairwise_best_auc.png",
        "family_median_heatmap.png",
        "pca_families.png",
        "pca_length_adjusted.png",
        "box_top_features.png",
        "box_residual_features.png",
        "box_uncorrelated_features.png",
        "confusion_ridge.png",
    ):
        print("wrote %s" % (output_dir / name))


if __name__ == "__main__":
    main()
