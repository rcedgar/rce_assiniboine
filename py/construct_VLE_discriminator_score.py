#!/usr/bin/env python3
"""Fit a closed-form VLE score and a precision-0.5 operating point.

S is a nonnegative weighted sum of robust z-scores. Each training fold
multiplies a feature by -1 when that direction has the better recall at
precision >= 0.5. Counts and lengths are log1p before the z-score. nr_bases
and exact duplicate columns stay out. After a term is accepted, a remaining
candidate whose family-weighted absolute Pearson correlation with it is at
least 0.90 is dropped. Forward selection adds a term only while
family-weighted recall at precision >= 0.5 rises. The number of terms is the
median length chosen by leave-one-VLE-family-out. The published formula is
fit on every row. Validation holds out one VLE family at a time.
"""

from __future__ import annotations

import math
import multiprocessing
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
_DEFAULT_OUTPUT = _ROOT / "construct_VLE_discriminator_score"
_TARGET_PRECISION = 0.5
_POLISH_EVALS = 120
_MIN_FAMILY_N = 20
_LENGTH_RHO_MAX = 0.30
_CORR_DROP = 0.90
_LENGTH_FEATURE = "nr_bases"
_MIX_ALPHAS = (0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85)
_STEM_FEATURE = "longest_composite_stem_total_nr_bases"
_COLOR_VLE = "#d55e00"
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
    "pOffRod",
    "pOnRod",
    "pSplitStem",
    "pLongestPairShare",
    "pJunction",
    "pDominantSpoke",
}
_WORKER: dict = {}


class _Point(NamedTuple):
    recall: float
    precision: float
    fpr: float
    threshold: float


class _Step(NamedTuple):
    feature: int
    weights: np.ndarray
    train_point: _Point


class _View(NamedTuple):
    z: np.ndarray
    is_vle: np.ndarray
    row_weight: np.ndarray
    median: np.ndarray
    iqr: np.ndarray
    active: np.ndarray
    sign: np.ndarray


class _Search(NamedTuple):
    steps: tuple[_Step, ...]
    median: np.ndarray
    iqr: np.ndarray
    active: np.ndarray
    rho: np.ndarray
    eligible: tuple[int, ...]
    drops: tuple[tuple[str, str, float], ...]
    sign: np.ndarray


class _Outer(NamedTuple):
    family: str
    k: int
    n_test: int
    n_vle: int
    holdout_called: int
    stem_called: int
    train_point: _Point
    weights: np.ndarray
    path_features: tuple[int, ...]
    inner_names: tuple[str, ...]
    inner_recalls: tuple[tuple[float, ...], ...]
    gain_rows: tuple[tuple[int, float, float, bool], ...]
    voting: bool
    auc: float
    stem_auc: float
    fpr: np.ndarray
    tpr: np.ndarray


class _Model(NamedTuple):
    weights: np.ndarray
    median: np.ndarray
    iqr: np.ndarray
    train_point: _Point
    n_eval: int
    k: int
    requested_k: int
    sign: np.ndarray


class _Candidate(NamedTuple):
    recall: float
    precision: float
    feature: int
    alpha: float
    weights: np.ndarray
    point: _Point
    peak_precision: float
    peak_recall: float


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _expect(condition: bool, message: str) -> None:
    if not condition:
        _die("ERROR: self-check failed: " + message)


def _fmt(value: float) -> str:
    if not np.isfinite(value):
        return "inf" if value > 0.0 else "-inf"
    return "%.6g" % value


def _canonical_feature_names() -> list[str]:
    py_dir = str(Path(__file__).resolve().parent)
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)
    from dotbracket_circular_features import extract_circular_rna_features

    names = sorted(extract_circular_rna_features("").keys())
    if not names:
        _die("ERROR: extract_circular_rna_features returned no features")
    return names


def _is_header(first_field: str) -> bool:
    return first_field in ("#label", "label")


def _transform_of(name: str) -> str:
    if name.startswith("pct_") or "_pct_" in name or name in _IDENTITY_NAMES:
        return "identity"
    return "log1p"


def _terms() -> list[tuple[str, str]]:
    return [(name, _transform_of(name)) for name in _canonical_feature_names()]


def _family_of(label: str, line_number: int) -> str:
    parts = label.split("_", 2)
    if len(parts) != 3 or parts[0] not in ("VLE", "DECOY") or not parts[1] or not parts[2]:
        _die(
            "ERROR: line %d label '%s' is not VLE_family_... or DECOY_family_..."
            % (line_number, label)
        )
    return parts[1]


def _class_is_vle(label: str, line_number: int) -> bool:
    prefix = label.split("_", 1)[0]
    if prefix == "VLE":
        return True
    if prefix == "DECOY":
        return False
    _die("ERROR: line %d label '%s' is not VLE or DECOY" % (line_number, label))


def _load(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[tuple[str, str]]]:
    """Return every feature, a VLE mask, the family of each row, and the term list."""
    terms = _terms()
    wanted = [name for name, _kind in terms]
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
            missing = [name for name in wanted if name not in feature_names]
            if missing:
                _die("ERROR: %s is missing %s" % (path, ", ".join(missing)))
            columns = [feature_names.index(name) + 1 for name in wanted]
            width = 1 + len(feature_names)
            raw_columns: list[list[float]] = [[] for _name in wanted]
            is_vle_list: list[bool] = []
            family_list: list[str] = []

            def consume(line_number: int, fields: list[str]) -> None:
                if len(fields) != width:
                    _die(
                        "ERROR: %s:%d has %d columns, expected %d"
                        % (path, line_number, len(fields), width)
                    )
                is_vle_list.append(_class_is_vle(fields[0], line_number))
                family_list.append(_family_of(fields[0], line_number))
                for slot, column in enumerate(columns):
                    try:
                        value = float(fields[column])
                    except ValueError:
                        _die(
                            "ERROR: %s:%d %s is not a number"
                            % (path, line_number, wanted[slot])
                        )
                    if not np.isfinite(value):
                        _die(
                            "ERROR: %s:%d %s is not finite"
                            % (path, line_number, wanted[slot])
                        )
                    raw_columns[slot].append(value)

            if pending_first is not None:
                consume(1, pending_first)
            next_report = 100000
            for line_number, line in enumerate(handle, start=2):
                text = line.rstrip("\r\n")
                if text == "":
                    continue
                consume(line_number, text.split("\t"))
                if len(is_vle_list) >= next_report:
                    print("loaded %d rows" % len(is_vle_list), file=sys.stderr, flush=True)
                    next_report += 100000
    except OSError as exc:
        _die("ERROR: cannot read '%s': %s" % (path, exc.strerror or exc))
    if not is_vle_list:
        _die("ERROR: no data rows in %s" % path)
    raw = np.column_stack([np.asarray(column, dtype=np.float64) for column in raw_columns])
    is_vle = np.asarray(is_vle_list, dtype=bool)
    families = np.asarray(family_list)
    if not bool(is_vle.any()) or not bool((~is_vle).any()):
        _die("ERROR: need both VLE and DECOY rows")
    if np.any(families == "avsunviroid"):
        _die("ERROR: %s contains avsunviroid" % path)
    for index, (name, kind) in enumerate(terms):
        if kind == "log1p" and np.any(raw[:, index] < 0.0):
            _die("ERROR: %s is negative" % name)
    return raw, is_vle, families, terms


def _transform(raw: np.ndarray, terms: list[tuple[str, str]]) -> np.ndarray:
    transformed = np.empty_like(raw)
    for index, (_name, kind) in enumerate(terms):
        if kind == "log1p":
            transformed[:, index] = np.log1p(raw[:, index])
        else:
            transformed[:, index] = raw[:, index]
    return transformed


def _fit_z(train_raw: np.ndarray, terms: list[tuple[str, str]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    transformed = _transform(train_raw, terms)
    median = np.median(transformed, axis=0)
    q25, q75 = np.quantile(transformed, [0.25, 0.75], axis=0)
    iqr = q75 - q25
    active = iqr > 0.0
    if not np.any(active):
        _die("ERROR: every training feature has IQR 0")
    return median, iqr, active


def _apply_z(
    raw: np.ndarray,
    median: np.ndarray,
    iqr: np.ndarray,
    active: np.ndarray,
    terms: list[tuple[str, str]],
) -> np.ndarray:
    transformed = _transform(raw, terms)
    scaled = np.zeros_like(transformed)
    scaled[:, active] = (transformed[:, active] - median[active]) / iqr[active]
    return scaled


def _signed_z(raw_rows: np.ndarray, view: _View) -> np.ndarray:
    """Robust z-scores with this training fold's feature signs applied."""
    return _apply_z(raw_rows, view.median, view.iqr, view.active, _WORKER["terms"]) * view.sign


def _orient_columns(
    z: np.ndarray,
    sign: np.ndarray,
    columns: list[int],
    is_vle: np.ndarray,
    row_weight: np.ndarray,
) -> None:
    """Flip a column when -z wins the training operating point. A tie stays positive."""
    weights = np.zeros(z.shape[1], dtype=np.float64)
    for feature in columns:
        column = z[:, feature]
        pos_point = _point_at_precision(column, is_vle, _TARGET_PRECISION, row_weight)
        neg_point = _point_at_precision(-column, is_vle, _TARGET_PRECISION, row_weight)
        pos_peak_precision, pos_peak_recall = _peak(column, is_vle, row_weight)
        neg_peak_precision, neg_peak_recall = _peak(-column, is_vle, row_weight)
        positive = _Candidate(
            pos_point.recall, pos_point.precision, feature, 1.0, weights, pos_point,
            pos_peak_precision, pos_peak_recall,
        )
        negative = _Candidate(
            neg_point.recall, neg_point.precision, feature, 1.0, weights, neg_point,
            neg_peak_precision, neg_peak_recall,
        )
        if _prefer(negative, positive):
            sign[feature] = -1.0
            z[:, feature] = -column


def _weighted_pearson_columns(
    column: np.ndarray, others: np.ndarray, row_weight: np.ndarray
) -> np.ndarray:
    """Family-weighted Pearson of one column against each column of others."""
    count = int(others.shape[1])
    corr = np.zeros(count, dtype=np.float64)
    total = float(row_weight.sum())
    if total <= 0.0 or count == 0:
        return corr
    centered_x = column - float(np.dot(row_weight, column) / total)
    var_x = float(np.dot(row_weight, centered_x * centered_x) / total)
    if var_x <= 0.0:
        return corr
    means = others.T @ row_weight / total
    centered = others - means
    weighted = centered * row_weight[:, None]
    cov = weighted.T @ centered_x / total
    var_y = (weighted * centered).sum(axis=0) / total
    ok = var_y > 0.0
    corr[ok] = cov[ok] / np.sqrt(var_x * var_y[ok])
    return corr


def _tie_ends(sorted_scores: np.ndarray) -> np.ndarray:
    """Last index of each tied run in a score array sorted in either direction."""
    change = np.empty(sorted_scores.size, dtype=bool)
    change[-1] = True
    change[:-1] = sorted_scores[:-1] != sorted_scores[1:]
    return np.flatnonzero(change)


def _average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="quicksort")
    sorted_values = values[order]
    change = np.empty(values.size, dtype=bool)
    change[0] = True
    change[1:] = sorted_values[1:] != sorted_values[:-1]
    starts = np.flatnonzero(change)
    ends = np.empty_like(starts)
    ends[:-1] = starts[1:] - 1
    ends[-1] = values.size - 1
    average = 0.5 * ((starts + 1) + (ends + 1))
    ranks = np.empty(values.size, dtype=np.float64)
    ranks[order] = np.repeat(average, ends - starts + 1)
    return ranks


def _ranked_curves(
    scores: np.ndarray, is_vle: np.ndarray, row_weight: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Tie-end weighted recall, precision, FPR, and the score of each group."""
    pos_mass = float(row_weight[is_vle].sum())
    neg_mass = float(row_weight[~is_vle].sum())
    if pos_mass <= 0.0 or neg_mass <= 0.0:
        _die("ERROR: a split is missing a class")
    order = np.argsort(-scores, kind="quicksort")
    sorted_scores = scores[order]
    ordered_pos = is_vle[order]
    ordered_weight = row_weight[order]
    true_positive = np.cumsum(np.where(ordered_pos, ordered_weight, 0.0))
    false_positive = np.cumsum(np.where(ordered_pos, 0.0, ordered_weight))
    ends = _tie_ends(sorted_scores)
    tp = true_positive[ends]
    fp = false_positive[ends]
    predicted = tp + fp
    recall = tp / pos_mass
    precision = np.zeros_like(tp)
    np.divide(tp, predicted, out=precision, where=predicted > 0.0)
    fpr = fp / neg_mass
    return recall, precision, fpr, sorted_scores[ends]


def _point_at_precision(
    scores: np.ndarray,
    is_vle: np.ndarray,
    target: float,
    row_weight: np.ndarray | None = None,
) -> _Point:
    """Largest weighted recall whose threshold still has precision >= target.

    Ties share one threshold. Among equal recalls, the higher threshold wins.
    """
    if row_weight is None:
        row_weight = np.ones(scores.size, dtype=np.float64)
    recall, precision, fpr, group_score = _ranked_curves(scores, is_vle, row_weight)
    ok = precision >= target
    if not np.any(ok):
        return _Point(0.0, 0.0, 0.0, np.inf)
    best_recall = float(recall[ok].max())
    chosen = int(np.flatnonzero(ok & (recall >= best_recall - 1e-12))[0])
    return _Point(
        float(recall[chosen]),
        float(precision[chosen]),
        float(fpr[chosen]),
        float(group_score[chosen]),
    )


def _auc(scores: np.ndarray, is_vle: np.ndarray) -> float:
    n_pos = int(is_vle.sum())
    n_neg = int(scores.size - n_pos)
    if n_pos == 0 or n_neg == 0 or np.all(scores == scores[0]):
        return 0.5
    order = np.argsort(scores, kind="quicksort")
    sorted_scores = scores[order]
    ends = _tie_ends(sorted_scores)
    starts = np.empty_like(ends)
    starts[0] = 0
    starts[1:] = ends[:-1] + 1
    average = 0.5 * ((starts + 1) + (ends + 1))
    ranks = np.empty(scores.size, dtype=np.float64)
    ranks[order] = np.repeat(average, ends - starts + 1)
    rank_sum = float(ranks[is_vle].sum())
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def _roc_points(scores: np.ndarray, is_vle: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unweighted ROC, one point per tied score, starting at the origin."""
    n_pos = int(is_vle.sum())
    n_neg = int(scores.size - n_pos)
    if n_pos == 0 or n_neg == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])
    order = np.argsort(-scores, kind="quicksort")
    sorted_scores = scores[order]
    ordered_pos = is_vle[order]
    true_positive = np.cumsum(ordered_pos)
    false_positive = np.cumsum(~ordered_pos)
    ends = _tie_ends(sorted_scores)
    fpr = np.empty(ends.size + 1, dtype=np.float64)
    tpr = np.empty(ends.size + 1, dtype=np.float64)
    fpr[0] = 0.0
    tpr[0] = 0.0
    fpr[1:] = false_positive[ends] / n_neg
    tpr[1:] = true_positive[ends] / n_pos
    if fpr.size <= 800:
        return fpr, tpr
    index = np.unique(np.linspace(0, fpr.size - 1, 800).astype(np.int64))
    return fpr[index], tpr[index]


def _holdout_ranking(
    view: _View, weights: np.ndarray, test_index: np.ndarray
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """AUC of one held-out family against every decoy, plus the stem column."""
    decoy_index = np.flatnonzero(~_WORKER["is_vle"])
    rows = np.concatenate([test_index, decoy_index])
    scored = _signed_z(_WORKER["raw"][rows], view)
    labels = np.zeros(rows.size, dtype=bool)
    labels[: test_index.size] = True
    formula = scored @ weights
    stem = scored[:, _WORKER["stem_index"]]
    return _auc(formula, labels), _auc(stem, labels), *_roc_points(formula, labels)


def _project(weights: np.ndarray) -> np.ndarray:
    clipped = np.maximum(weights, 0.0)
    total = float(clipped.sum())
    if total <= 0.0:
        return np.full(weights.size, 1.0 / weights.size)
    return clipped / total


def _project_support(weights: np.ndarray, support: np.ndarray) -> np.ndarray:
    clipped = np.zeros(weights.size, dtype=np.float64)
    clipped[support] = np.maximum(weights[support], 0.0)
    total = float(clipped[support].sum())
    if total <= 0.0:
        clipped[support] = 1.0 / support.size
        return clipped
    clipped[support] /= total
    return clipped


def _vle_family_names(families: np.ndarray, is_vle: np.ndarray) -> list[str]:
    names = []
    for name in np.unique(families[is_vle]):
        member = families == name
        if bool(np.all(is_vle[member])):
            names.append(str(name))
    return sorted(names)


def _voting_names(families: np.ndarray, is_vle: np.ndarray) -> list[str]:
    names = []
    for name in _vle_family_names(families, is_vle):
        count = int(np.count_nonzero((families == name) & is_vle))
        if count >= _MIN_FAMILY_N:
            names.append(name)
    return names


def _row_weights(families: np.ndarray, is_vle: np.ndarray) -> np.ndarray:
    """Equal total weight per large VLE family and per DECOY family.

    Families under the voting cutoff stay at weight 1 and sit outside the VLE
    split. Each class's weights sum to that class's row count.
    """
    weights = np.ones(is_vle.size, dtype=np.float64)
    voting: list[tuple[np.ndarray, int]] = []
    small_n = 0
    for name in _vle_family_names(families, is_vle):
        index = np.flatnonzero((families == name) & is_vle)
        count = int(index.size)
        if count >= _MIN_FAMILY_N:
            voting.append((index, count))
        else:
            small_n += count
    if voting:
        mass = float(int(is_vle.sum()) - small_n)
        share = mass / len(voting)
        for index, count in voting:
            weights[index] = share / count
    decoy_names = sorted({str(name) for name in families[~is_vle]})
    if decoy_names:
        decoy_share = float(int((~is_vle).sum())) / len(decoy_names)
        for name in decoy_names:
            index = np.flatnonzero((families == name) & ~is_vle)
            count = int(index.size)
            if count == 0:
                continue
            weights[index] = decoy_share / count
    return weights


def _column_rho(train_raw: np.ndarray, length_index: int) -> np.ndarray:
    reference = _average_ranks(train_raw[:, length_index])
    centered_ref = reference - reference.mean()
    ref_norm = float(np.sqrt(np.dot(centered_ref, centered_ref)))
    rho = np.zeros(train_raw.shape[1], dtype=np.float64)
    for index in range(train_raw.shape[1]):
        ranks = _average_ranks(train_raw[:, index])
        centered = ranks - ranks.mean()
        norm = float(np.sqrt(np.dot(centered, centered)))
        if ref_norm == 0.0 or norm == 0.0:
            rho[index] = 0.0
        else:
            rho[index] = float(np.dot(centered_ref, centered) / (ref_norm * norm))
    return rho


def _eligible_columns(
    train_raw: np.ndarray,
    terms: list[tuple[str, str]],
    active: np.ndarray,
    rho: np.ndarray,
) -> tuple[list[int], list[tuple[str, str, float]]]:
    kept: list[int] = []
    drops: list[tuple[str, str, float]] = []
    candidates: list[int] = []
    for index, (name, _kind) in enumerate(terms):
        value = float(rho[index])
        if name == _LENGTH_FEATURE:
            drops.append((name, "length", value))
            continue
        if not bool(active[index]):
            drops.append((name, "iqr0", value))
            continue
        candidates.append(index)
    for index in candidates:
        name = terms[index][0]
        duplicate_of = None
        for previous in kept:
            if np.array_equal(train_raw[:, index], train_raw[:, previous]):
                duplicate_of = terms[previous][0]
                break
        if duplicate_of is not None:
            drops.append((name, "duplicate of %s" % duplicate_of, float(rho[index])))
            continue
        kept.append(index)
    return kept, drops


def _peak(scores: np.ndarray, is_vle: np.ndarray, row_weight: np.ndarray) -> tuple[float, float]:
    """Highest precision on the curve, and the recall at that precision."""
    recall, precision, _fpr, _score = _ranked_curves(scores, is_vle, row_weight)
    chosen = int(precision.argmax())
    return float(precision[chosen]), float(recall[chosen])


def _prefer(cand: _Candidate, best: _Candidate | None) -> bool:
    if best is None:
        return True
    if cand.recall > best.recall + 1e-12:
        return True
    if cand.recall < best.recall - 1e-12:
        return False
    if cand.precision > best.precision + 1e-12:
        return True
    if cand.precision < best.precision - 1e-12:
        return False
    if cand.peak_precision > best.peak_precision + 1e-12:
        return True
    if cand.peak_precision < best.peak_precision - 1e-12:
        return False
    if cand.peak_recall > best.peak_recall + 1e-12:
        return True
    if cand.peak_recall < best.peak_recall - 1e-12:
        return False
    if cand.feature != best.feature:
        return cand.feature < best.feature
    return cand.alpha < best.alpha


def _mix_weights(current: np.ndarray, feature: int, alpha: float) -> np.ndarray:
    weights = current * (1.0 - alpha)
    weights[feature] = current[feature] * (1.0 - alpha) + alpha
    return weights


def _evaluate_mix(
    base_score: np.ndarray,
    column: np.ndarray,
    current: np.ndarray,
    feature: int,
    alpha: float,
    is_vle: np.ndarray,
    row_weight: np.ndarray,
) -> _Candidate:
    if float(np.abs(current).sum()) == 0.0:
        score = column
        weights = np.zeros(current.size, dtype=np.float64)
        weights[feature] = 1.0
        alpha = 1.0
    else:
        score = (1.0 - alpha) * base_score + alpha * column
        weights = _mix_weights(current, feature, alpha)
    point = _point_at_precision(score, is_vle, _TARGET_PRECISION, row_weight)
    peak_precision, peak_recall = _peak(score, is_vle, row_weight)
    return _Candidate(
        point.recall,
        point.precision,
        feature,
        alpha,
        weights,
        point,
        peak_precision,
        peak_recall,
    )


def _forward_path(
    z: np.ndarray,
    is_vle: np.ndarray,
    row_weight: np.ndarray,
    eligible: list[int],
    terms: list[tuple[str, str]],
    label: str,
) -> tuple[tuple[_Step, ...], tuple[tuple[str, str, float], ...]]:
    steps: list[_Step] = []
    corr_drops: list[tuple[str, str, float]] = []
    current = np.zeros(z.shape[1], dtype=np.float64)
    current_recall = -1.0
    remaining = list(eligible)
    while remaining:
        best: _Candidate | None = None
        base_score = z @ current
        opening = float(np.abs(current).sum()) == 0.0
        for feature in remaining:
            column = z[:, feature]
            alphas = (1.0,) if opening else _MIX_ALPHAS
            for alpha in alphas:
                cand = _evaluate_mix(
                    base_score, column, current, feature, alpha, is_vle, row_weight
                )
                if opening and label == "all-data":
                    print(
                        "  %s  screen %s  recall %.4f  precision %.4f  peak_precision %.4f"
                        % (
                            label,
                            terms[feature][0],
                            cand.recall,
                            cand.precision,
                            cand.peak_precision,
                        ),
                        file=sys.stderr,
                        flush=True,
                    )
                if _prefer(cand, best):
                    best = cand
        if best is None or best.recall <= current_recall + 1e-12:
            break
        steps.append(_Step(best.feature, best.weights.copy(), best.point))
        current = best.weights
        current_recall = best.recall
        chosen = best.feature
        remaining.remove(chosen)
        if remaining:
            indexes = np.array(remaining, dtype=np.int64)
            corrs = _weighted_pearson_columns(z[:, chosen], z[:, indexes], row_weight)
            kept: list[int] = []
            chosen_name = terms[chosen][0]
            for feature, corr in zip(remaining, corrs):
                value = float(corr)
                if abs(value) >= _CORR_DROP:
                    name = terms[feature][0]
                    corr_drops.append((name, "correlated with %s" % chosen_name, value))
                    print(
                        "  %s  drop %s  correlated with %s  pearson %.3f"
                        % (label, name, chosen_name, value),
                        file=sys.stderr,
                        flush=True,
                    )
                else:
                    kept.append(feature)
            remaining = kept
        print(
            "  %s  term %d  %s  alpha %.2f  recall %.4f"
            % (label, len(steps), terms[best.feature][0], best.alpha, best.recall),
            file=sys.stderr,
            flush=True,
        )
    if not steps:
        _die("ERROR: %s has no eligible feature" % label)
    return tuple(steps), tuple(corr_drops)


def _refit_order(
    z: np.ndarray,
    is_vle: np.ndarray,
    row_weight: np.ndarray,
    order: list[int] | tuple[int, ...],
) -> list[tuple[np.ndarray, _Point]]:
    current = np.zeros(z.shape[1], dtype=np.float64)
    fitted: list[tuple[np.ndarray, _Point]] = []
    for feature in order:
        best: _Candidate | None = None
        base_score = z @ current
        column = z[:, feature]
        opening = float(np.abs(current).sum()) == 0.0
        alphas = (1.0,) if opening else _MIX_ALPHAS
        for alpha in alphas:
            cand = _evaluate_mix(
                base_score, column, current, feature, alpha, is_vle, row_weight
            )
            if _prefer(cand, best):
                best = cand
        if best is None:
            _die("ERROR: a fixed feature order produced no weights")
        current = best.weights
        fitted.append((current.copy(), best.point))
    return fitted


class _Budget:
    def __init__(self, limit: int, row_weight: np.ndarray):
        self.limit = limit
        self.row_weight = row_weight
        self.n_eval = 0
        self.cache: dict[tuple[float, ...], float] = {}
        self.exhausted = False

    def recall(self, weights: np.ndarray, features: np.ndarray, is_vle: np.ndarray) -> float | None:
        key = tuple(np.round(weights, 8))
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        if self.n_eval >= self.limit:
            self.exhausted = True
            return None
        self.n_eval += 1
        point = _point_at_precision(features @ weights, is_vle, _TARGET_PRECISION, self.row_weight)
        self.cache[key] = point.recall
        return point.recall


def _exploratory(
    origin: np.ndarray,
    origin_value: float,
    delta: float,
    budget: _Budget,
    features: np.ndarray,
    is_vle: np.ndarray,
    support: np.ndarray,
) -> tuple[np.ndarray, float]:
    current = origin.copy()
    current_value = origin_value
    for coordinate in support:
        best = current
        best_value = current_value
        for sign in (1.0, -1.0):
            candidate = current.copy()
            candidate[coordinate] += sign * delta
            candidate = _project_support(candidate, support)
            value = budget.recall(candidate, features, is_vle)
            if value is not None and value > best_value + 1e-12:
                best = candidate
                best_value = value
        current = best
        current_value = best_value
    return current, current_value


def _hooke_jeeves(
    start: np.ndarray,
    budget: _Budget,
    features: np.ndarray,
    is_vle: np.ndarray,
    support: np.ndarray,
) -> tuple[np.ndarray, float]:
    base = _project_support(start, support)
    base_value = budget.recall(base, features, is_vle)
    if base_value is None:
        _die("ERROR: Hooke-Jeeves has no evaluations left at the start")
    delta = 0.25
    while delta >= 0.01 and not budget.exhausted:
        trial, trial_value = _exploratory(
            base, base_value, delta, budget, features, is_vle, support
        )
        if trial_value <= base_value + 1e-12:
            delta *= 0.5
            continue
        previous = base
        base = trial
        base_value = trial_value
        jumps = 0
        while jumps < 25 and not budget.exhausted:
            jumps += 1
            pattern = _project_support(base + (base - previous), support)
            pattern_value = budget.recall(pattern, features, is_vle)
            if pattern_value is None:
                break
            anchor = pattern_value if pattern_value > base_value else base_value
            refined, refined_value = _exploratory(
                pattern, anchor, delta, budget, features, is_vle, support
            )
            if refined_value > base_value + 1e-12:
                previous = base
                base = refined
                base_value = refined_value
            else:
                break
    return base, base_value


def _polish(
    z: np.ndarray,
    is_vle: np.ndarray,
    row_weight: np.ndarray,
    start: np.ndarray,
    support: np.ndarray,
) -> tuple[np.ndarray, _Point, int]:
    budget = _Budget(_POLISH_EVALS, row_weight)
    weights, _recall = _hooke_jeeves(start, budget, z, is_vle, support)
    point = _point_at_precision(z @ weights, is_vle, _TARGET_PRECISION, row_weight)
    return weights, point, budget.n_eval


def _positive_calls(scores: np.ndarray, is_vle: np.ndarray, threshold: float) -> tuple[int, int]:
    n_pos = int(is_vle.sum())
    if n_pos == 0 or not np.isfinite(threshold):
        return 0, n_pos
    called = int(np.count_nonzero((scores >= threshold) & is_vle))
    return called, n_pos


def _choose_k(
    recall_rows: list[np.ndarray],
) -> tuple[int, tuple[tuple[int, float, float, bool], ...]]:
    """Keep adding a term while the mean family gain exceeds its standard error."""
    if not recall_rows:
        _die("ERROR: a complexity path has no prefix")
    kept = 1
    notes: list[tuple[int, float, float, bool]] = []
    for k in range(2, len(recall_rows) + 1):
        gain = recall_rows[k - 1] - recall_rows[k - 2]
        count = int(gain.size)
        mean = float(gain.mean()) if count else 0.0
        if count < 2:
            notes.append((k, mean, math.inf, False))
            break
        sd = float(gain.std(ddof=1))
        se = sd / math.sqrt(count)
        accepted = mean > se
        notes.append((k, mean, se, accepted))
        if not accepted:
            break
        kept = k
    return kept, tuple(notes)


def _median_k(values: list[int]) -> int:
    if not values:
        _die("ERROR: no voting family chose a formula length")
    ordered = sorted(values)
    count = len(ordered)
    if count % 2 == 1:
        return int(ordered[count // 2])
    average = 0.5 * (ordered[count // 2 - 1] + ordered[count // 2])
    return int(math.floor(average))


def _mask_without(names: list[str]) -> np.ndarray:
    mask = np.ones(_WORKER["n"], dtype=bool)
    rows = _WORKER["family_rows"]
    for name in names:
        mask[rows[name]] = False
    return mask


def _training_view(mask: np.ndarray) -> _View:
    raw = _WORKER["raw"]
    terms = _WORKER["terms"]
    train_raw = raw[mask]
    median, iqr, active = _fit_z(train_raw, terms)
    z = _apply_z(train_raw, median, iqr, active, terms)
    is_vle = _WORKER["is_vle"][mask]
    row_weight = _row_weights(_WORKER["families"][mask], is_vle)
    sign = np.ones(z.shape[1], dtype=np.float64)
    return _View(z, is_vle, row_weight, median, iqr, active, sign)


def _search_on_view(view: _View, mask: np.ndarray, label: str) -> _Search:
    terms = _WORKER["terms"]
    train_raw = _WORKER["raw"][mask]
    rho = _column_rho(train_raw, _WORKER["length_index"])
    eligible, drops = _eligible_columns(train_raw, terms, view.active, rho)
    if not eligible:
        _die("ERROR: %s has no eligible feature" % label)
    print(
        "%s  eligible %d  %s"
        % (label, len(eligible), ", ".join(terms[index][0] for index in eligible)),
        file=sys.stderr,
        flush=True,
    )
    if label == "all-data":
        for name, reason, value in drops:
            print(
                "  %s  drop %s  %s  rho %.3f" % (label, name, reason, value),
                file=sys.stderr,
                flush=True,
            )
    columns = list(eligible)
    stem_index = _WORKER["stem_index"]
    if stem_index not in columns:
        columns.append(stem_index)
    _orient_columns(view.z, view.sign, columns, view.is_vle, view.row_weight)
    flipped = [terms[index][0] for index in columns if view.sign[index] < 0.0]
    print(
        "%s  flipped %d  %s" % (label, len(flipped), ", ".join(flipped) if flipped else "-"),
        file=sys.stderr,
        flush=True,
    )
    steps, corr_drops = _forward_path(view.z, view.is_vle, view.row_weight, eligible, terms, label)
    return _Search(
        steps,
        view.median,
        view.iqr,
        view.active,
        rho,
        tuple(eligible),
        tuple(drops) + corr_drops,
        view.sign.copy(),
    )


def _stem_called(view: _View, train_mask: np.ndarray, test_index: np.ndarray) -> int:
    stem_index = _WORKER["stem_index"]
    raw = _WORKER["raw"]
    train_z = _signed_z(raw[train_mask], view)
    point = _point_at_precision(
        train_z[:, stem_index], view.is_vle, _TARGET_PRECISION, view.row_weight
    )
    test_z = _signed_z(raw[test_index], view)
    called, _n_pos = _positive_calls(
        test_z[:, stem_index], _WORKER["is_vle"][test_index], point.threshold
    )
    return called


def _score_family(
    view: _View, weights: np.ndarray, threshold: float, test_index: np.ndarray
) -> tuple[int, int]:
    test_z = _signed_z(_WORKER["raw"][test_index], view)
    return _positive_calls(test_z @ weights, _WORKER["is_vle"][test_index], threshold)


def _orient_order(view: _View, order: tuple[int, ...]) -> None:
    """Orient the refit columns and the stem on this fold's own rows."""
    columns = list(order)
    stem_index = _WORKER["stem_index"]
    if stem_index not in columns:
        columns.append(stem_index)
    _orient_columns(view.z, view.sign, columns, view.is_vle, view.row_weight)


def _prefix_recalls(train_mask: np.ndarray, test_index: np.ndarray, order: tuple[int, ...]) -> list[float]:
    view = _training_view(train_mask)
    _orient_order(view, order)
    fitted = _refit_order(view.z, view.is_vle, view.row_weight, order)
    recalls = []
    for weights, point in fitted:
        called, n_pos = _score_family(view, weights, point.threshold, test_index)
        recalls.append((called / n_pos) if n_pos else 0.0)
    return recalls


def _fit_outer(family: str) -> _Outer:
    started = time.perf_counter()
    train_mask = _mask_without([family])
    test_index = _WORKER["family_rows"][family]
    if int((~_WORKER["is_vle"][test_index]).sum()) != 0:
        _die("ERROR: %s is not all VLE" % family)
    view = _training_view(train_mask)
    search = _search_on_view(view, train_mask, family)
    order = tuple(step.feature for step in search.steps)
    inner_names = tuple(name for name in _WORKER["voting"] if name != family)
    columns: list[list[float]] = []
    for name in inner_names:
        inner_train = train_mask.copy()
        inner_train[_WORKER["family_rows"][name]] = False
        columns.append(_prefix_recalls(inner_train, _WORKER["family_rows"][name], order))
    recall_rows = []
    inner_recalls: list[tuple[float, ...]] = []
    for k in range(len(order)):
        row = np.array([columns[item][k] for item in range(len(inner_names))], dtype=np.float64)
        recall_rows.append(row)
        inner_recalls.append(tuple(float(value) for value in row))
    chosen, notes = _choose_k(recall_rows)
    step = search.steps[chosen - 1]
    called, n_pos = _score_family(view, step.weights, step.train_point.threshold, test_index)
    stem_called = _stem_called(view, train_mask, test_index)
    auc, stem_auc, fpr, tpr = _holdout_ranking(view, step.weights, test_index)
    print(
        "family %s done  k %d  recall %.4f  auc %.4f  in %.1fs"
        % (family, chosen, (called / n_pos) if n_pos else 0.0, auc, time.perf_counter() - started),
        file=sys.stderr,
        flush=True,
    )
    return _Outer(
        family=family,
        k=chosen,
        n_test=int(test_index.size),
        n_vle=n_pos,
        holdout_called=called,
        stem_called=stem_called,
        train_point=step.train_point,
        weights=step.weights,
        path_features=order,
        inner_names=inner_names,
        inner_recalls=tuple(inner_recalls),
        gain_rows=notes,
        voting=True,
        auc=auc,
        stem_auc=stem_auc,
        fpr=fpr,
        tpr=tpr,
    )


def _fit_published_path() -> _Search:
    mask = np.ones(_WORKER["n"], dtype=bool)
    view = _training_view(mask)
    return _search_on_view(view, mask, "all-data")


def _holdout_fixed(family: str, feature_order: tuple[int, ...], k: int) -> _Outer:
    order = feature_order[:k]
    if not order:
        _die("ERROR: published formula has no terms for %s" % family)
    train_mask = _mask_without([family])
    test_index = _WORKER["family_rows"][family]
    view = _training_view(train_mask)
    _orient_order(view, order)
    fitted = _refit_order(view.z, view.is_vle, view.row_weight, order)
    weights, point = fitted[-1]
    called, n_pos = _score_family(view, weights, point.threshold, test_index)
    auc, stem_auc, fpr, tpr = _holdout_ranking(view, weights, test_index)
    return _Outer(
        family=family,
        k=k,
        n_test=int(test_index.size),
        n_vle=n_pos,
        holdout_called=called,
        stem_called=_stem_called(view, train_mask, test_index),
        train_point=point,
        weights=weights,
        path_features=tuple(order),
        inner_names=(),
        inner_recalls=(),
        gain_rows=(),
        voting=False,
        auc=auc,
        stem_auc=stem_auc,
        fpr=fpr,
        tpr=tpr,
    )


def _run_job(spec: tuple[str, str]) -> tuple[str, object]:
    kind, family = spec
    if kind == "path":
        return ("path", _fit_published_path())
    return ("outer", _fit_outer(family))


def _init_worker(
    raw: np.ndarray,
    is_vle: np.ndarray,
    families: np.ndarray,
    terms: list[tuple[str, str]],
    stem_index: int,
    length_index: int,
    voting: list[str],
    family_rows: dict[str, np.ndarray],
) -> None:
    _WORKER.clear()
    _WORKER["raw"] = raw
    _WORKER["is_vle"] = is_vle
    _WORKER["families"] = families
    _WORKER["terms"] = terms
    _WORKER["stem_index"] = stem_index
    _WORKER["length_index"] = length_index
    _WORKER["voting"] = list(voting)
    _WORKER["family_rows"] = family_rows
    _WORKER["n"] = int(raw.shape[0])


def _self_check() -> None:
    scores = np.array([5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    is_vle = np.array([True, True, False, True, False, False])
    point = _point_at_precision(scores, is_vle, 0.5)
    _expect(abs(point.recall - 1.0) < 1e-12, "full recall")
    _expect(abs(point.threshold - 2.0) < 1e-12, "threshold includes the last VLE")
    _expect(abs(point.precision - 0.75) < 1e-12, "precision at that threshold")

    scores = np.array([5.0, 4.0, 1.0, 0.0, -1.0])
    is_vle = np.array([True, False, False, False, True])
    point = _point_at_precision(scores, is_vle, 0.5)
    _expect(abs(point.recall - 0.5) < 1e-12, "recall stops when precision falls")
    _expect(abs(point.threshold - 5.0) < 1e-12, "higher threshold wins a recall tie")

    scores = np.array([2.0, 2.0, 1.0])
    is_vle = np.array([True, False, True])
    point = _point_at_precision(scores, is_vle, 0.5)
    _expect(abs(point.recall - 1.0) < 1e-12, "tie group is accepted together")
    _expect(abs(point.threshold - 1.0) < 1e-12, "later group reaches full recall")
    _expect(abs(_auc(np.array([1.0, 0.0]), np.array([True, False])) - 1.0) < 1e-12, "perfect AUC")
    _expect(abs(_auc(np.array([1.0, 1.0]), np.array([True, False])) - 0.5) < 1e-12, "tied AUC")
    perfect_fpr, perfect_tpr = _roc_points(np.array([1.0, 0.0]), np.array([True, False]))
    _expect(perfect_fpr[0] == 0.0 and perfect_tpr[0] == 0.0, "ROC starts at the origin")
    _expect(abs(perfect_fpr[1]) < 1e-12 and abs(perfect_tpr[1] - 1.0) < 1e-12, "a perfect score is recalled at zero FPR")
    _expect(abs(perfect_fpr[-1] - 1.0) < 1e-12 and abs(perfect_tpr[-1] - 1.0) < 1e-12, "ROC ends at (1, 1)")
    projected = _project(np.array([1.0, -2.0, 1.0, 0.0]))
    _expect(abs(float(projected.sum()) - 1.0) < 1e-12, "simplex sums to 1")
    _expect(float(projected.min()) >= -1e-12, "simplex is nonnegative")
    support = np.array([0, 2])
    supported = _project_support(np.array([1.0, 5.0, -2.0, 4.0]), support)
    _expect(abs(float(supported.sum()) - 1.0) < 1e-12, "support sums to 1")
    _expect(abs(float(supported[1]) + float(supported[3])) < 1e-12, "support leaves the other coordinates at 0")
    _expect(_transform_of("pct_paired_bases") == "identity", "percentage stays raw")
    _expect(_transform_of("pRod") == "identity", "pRod stays raw")
    _expect(_transform_of("pOffRod") == "identity", "pOffRod stays raw")
    _expect(_transform_of("pSplitStem") == "identity", "pSplitStem stays raw")
    _expect(_transform_of("pJunction") == "identity", "pJunction stays raw")
    _expect(_transform_of("pOnRod") == "identity", "pOnRod stays raw")
    _expect(_transform_of("pLongestPairShare") == "identity", "pLongestPairShare stays raw")
    _expect(_transform_of("pDominantSpoke") == "identity", "pDominantSpoke stays raw")
    _expect(_transform_of("nBranch") == "log1p", "branch count is log1p")
    _expect(_transform_of(_STEM_FEATURE) == "log1p", "stem count is log1p")
    _expect(_family_of("VLE_ambivirus_rest", 1) == "ambivirus", "family is the second field")
    _expect(_family_of("DECOY_ambivirus_rest", 1) == "ambivirus", "decoy family shares the name")

    heavy = np.array([6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    heavy_class = np.array([True, False, False, False, False, False, True])
    heavy_weight = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0])
    heavy_point = _point_at_precision(heavy, heavy_class, 0.5, heavy_weight)
    _expect(abs(heavy_point.recall - 0.25) < 1e-12, "family-style weights change recall")
    _expect(abs(heavy_point.threshold - 6.0) < 1e-12, "weighted threshold stays on the light VLE")
    blocked_scores = np.array([2.0, 1.0, 0.0])
    blocked_class = np.array([True, False, True])
    blocked_weight = np.array([1.0, 3.0, 1.0])
    blocked_point = _point_at_precision(blocked_scores, blocked_class, 0.5, blocked_weight)
    _expect(abs(blocked_point.recall - 0.5) < 1e-12, "a heavy decoy blocks the lower VLE")
    _expect(abs(blocked_point.threshold - 2.0) < 1e-12, "threshold stays above the heavy decoy")

    families = np.array(["big"] * 100 + ["small"] * 40 + ["tiny"] * 5 + ["wide"] * 4 + ["narrow"])
    class_mask = np.array([True] * 145 + [False] * 5)
    row_weight = _row_weights(families, class_mask)
    big_total = float(row_weight[:100].sum())
    small_total = float(row_weight[100:140].sum())
    _expect(abs(big_total - 70.0) < 1e-8, "voting families share the VLE mass equally")
    _expect(abs(small_total - big_total) < 1e-8, "the smaller voting family has the same total")
    _expect(abs(float(row_weight[140:145].sum()) - 5.0) < 1e-12, "family under 20 stays at weight 1")
    _expect(abs(float(row_weight[class_mask].sum()) - 145.0) < 1e-8, "VLE weights sum to the VLE count")
    _expect(abs(float(row_weight[145:149].sum()) - 2.5) < 1e-8, "decoy families share the decoy mass equally")
    _expect(abs(float(row_weight[149]) - 2.5) < 1e-8, "a one-row decoy family has the same total")
    _expect(abs(float(row_weight[~class_mask].sum()) - 5.0) < 1e-8, "DECOY weights sum to the DECOY count")

    mixed_families = np.array(["a", "a", "b", "c"])
    mixed_class = np.array([True, True, False, True])
    _expect(_vle_family_names(mixed_families, mixed_class) == ["a", "c"], "hold-out families are all VLE")
    _expect(_voting_names(mixed_families, mixed_class) == [], "families under 20 do not vote")

    rho_raw = np.column_stack(
        [
            np.array([1.0, 2.0, 3.0, 4.0]),
            np.array([1.0, 2.0, 3.0, 4.0]),
            np.array([4.0, 3.0, 2.0, 1.0]),
            np.array([1.0, 0.0, 1.0, 0.0]),
        ]
    )
    rho = _column_rho(rho_raw, 0)
    _expect(abs(rho[0] - 1.0) < 1e-8, "length correlates with itself")
    _expect(abs(rho[1] - 1.0) < 1e-8, "a copy correlates with length")
    _expect(abs(rho[2] + 1.0) < 1e-8, "a reversal correlates at -1")

    terms = [(_LENGTH_FEATURE, "log1p"), ("high", "identity"), ("keep", "identity"), ("copy", "identity")]
    table = np.column_stack(
        [
            np.array([1.0, 2.0, 3.0, 4.0]),
            np.array([0.0, 1.0, 0.0, 1.0]),
            np.array([1.0, 0.0, 0.0, 1.0]),
            np.array([1.0, 0.0, 0.0, 1.0]),
        ]
    )
    injected = np.array([1.0, _LENGTH_RHO_MAX, 0.1, 0.1])
    active = np.array([True, True, True, False])
    kept, drops = _eligible_columns(table, terms, active, injected)
    _expect(kept == [1, 2], "a length correlate with positive IQR stays eligible")
    reasons = {name: reason for name, reason, _value in drops}
    _expect(reasons[_LENGTH_FEATURE] == "length", "nr_bases is a length feature")
    _expect("high" not in reasons, "rho at the cutoff stays eligible")
    _expect(reasons["copy"] == "iqr0", "an inactive column is dropped before the duplicate check")

    active = np.array([True, True, True, True])
    kept, drops = _eligible_columns(table, terms, active, injected)
    _expect(kept == [1, 2], "nr_bases is dropped and the other active columns stay")
    reasons = {name: reason for name, reason, _value in drops}
    _expect(reasons["copy"] == "duplicate of keep", "an exact copy is dropped")

    recall_rows = [
        np.array([0.5, 0.5, 0.5, 0.5]),
        np.array([0.8, 0.6, 0.6, 0.6]),
        np.array([0.82, 0.6, 0.6, 0.6]),
    ]
    chosen, notes = _choose_k(recall_rows)
    _expect(chosen == 2, "one-SE keeps the gain that clears the family scatter")
    _expect(notes[0][0] == 2 and notes[0][3], "the second term is accepted")
    _expect(notes[1][0] == 3 and not notes[1][3], "the third term is rejected")
    _expect(_median_k([1, 2, 2, 3, 3]) == 2, "odd median")
    _expect(_median_k([1, 2]) == 1, "even median rounds down")
    zeros = np.zeros(1, dtype=np.float64)
    empty = _Point(0.0, 0.0, 0.0, np.inf)
    weaker = _Candidate(0.0, 0.0, 0, 1.0, zeros, empty, 0.10, 0.20)
    stronger = _Candidate(0.0, 0.0, 1, 1.0, zeros, empty, 0.40, 0.50)
    _expect(_prefer(stronger, weaker), "higher peak precision wins a recall tie")
    _expect(not _prefer(weaker, stronger), "lower peak precision loses a recall tie")

    decoy_high = np.array([0.0, 1.0, 5.0, 4.0, 3.0])
    decoy_labels = np.array([True, True, False, False, False])
    decoy_weight = np.ones(decoy_high.size)
    raw_point = _point_at_precision(decoy_high, decoy_labels, _TARGET_PRECISION, decoy_weight)
    _expect(abs(raw_point.recall) < 1e-12, "the raw decoy-higher column has no recall")
    oriented = decoy_high.reshape(-1, 1).copy()
    oriented_sign = np.ones(1, dtype=np.float64)
    _orient_columns(oriented, oriented_sign, [0], decoy_labels, decoy_weight)
    _expect(oriented_sign[0] < 0.0, "a decoy-higher column is stored as sign -1")
    _expect(bool(decoy_labels[int(np.argmax(oriented[:, 0]))]), "the oriented column ranks a VLE first")
    flipped_point = _point_at_precision(oriented[:, 0], decoy_labels, _TARGET_PRECISION, decoy_weight)
    _expect(abs(flipped_point.recall - 1.0) < 1e-12, "the oriented column recalls every VLE")

    tied = np.array([[1.0], [1.0], [-1.0], [-1.0]])
    tied_labels = np.array([True, False, True, False])
    tied_sign = np.ones(1, dtype=np.float64)
    _orient_columns(tied, tied_sign, [0], tied_labels, np.ones(4))
    _expect(tied_sign[0] > 0.0, "equal recall on both directions keeps +1")
    _expect(abs(float(tied[0, 0]) - 1.0) < 1e-12, "a tied column is not flipped")

    base = np.array([3.0, 3.0, 3.0, 3.0, 0.0, 0.0, 0.0, 0.0])
    ortho = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0])

    def _scale(target: float) -> float:
        return math.sqrt(2.25 * (1.0 / (target * target) - 1.0))

    partner = base + _scale(0.95) * ortho
    distant = base + _scale(0.50) * ortho
    paired = np.column_stack([base, partner, distant])
    pair_weight = np.ones(base.size)
    measured = _weighted_pearson_columns(paired[:, 0], paired[:, 1:], pair_weight)
    _expect(abs(float(measured[0]) - 0.95) < 1e-8, "the near copy correlates at 0.95")
    _expect(abs(float(measured[1]) - 0.50) < 1e-8, "the distant column correlates at 0.50")
    pair_labels = np.array([True, True, True, True, False, False, False, False])
    pair_terms = [("base", "identity"), ("close", "identity"), ("mid", "identity")]
    pair_steps, pair_drops = _forward_path(
        paired, pair_labels, pair_weight, [0, 1, 2], pair_terms, "self-check"
    )
    _expect(len(pair_steps) == 1 and pair_steps[0].feature == 0, "the separable column is selected")
    _expect(len(pair_drops) == 1 and pair_drops[0][0] == "close", "a correlation of 0.95 drops the partner")
    _expect(pair_drops[0][1] == "correlated with base", "the drop names the accepted term")
    _expect(all(name != "mid" for name, _reason, _value in pair_drops), "a correlation of 0.50 does not")


def _rate(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0.0
    return numerator / denominator


def _feature_names(terms: list[tuple[str, str]], indexes: tuple[int, ...]) -> str:
    if not indexes:
        return ""
    return "\t".join(terms[index][0] for index in indexes)


def _write_parameters(path: Path, terms: list[tuple[str, str]], model: _Model) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("feature\ttransform\ttrain_median\ttrain_iqr\tweight\tsign\n")
        for index, (name, transform) in enumerate(terms):
            handle.write(
                "%s\t%s\t%s\t%s\t%s\t%+d\n"
                % (
                    name,
                    transform,
                    _fmt(float(model.median[index])),
                    _fmt(float(model.iqr[index])),
                    _fmt(float(model.weights[index])),
                    int(model.sign[index]),
                )
            )
        handle.write("threshold\t\t\t\t%s\t\n" % _fmt(model.train_point.threshold))


def _write_holdouts(path: Path, terms: list[tuple[str, str]], holdouts: list[_Outer]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        weight_header = "\t".join(name for name, _kind in terms)
        handle.write(
            "family\tn_test\tn_vle\tk\t"
            "train_weighted_recall\ttrain_weighted_precision\ttrain_threshold\t"
            "holdout_recall\tstem_holdout_recall\t%s\n" % weight_header
        )
        for holdout in holdouts:
            weights = "\t".join(_fmt(float(value)) for value in holdout.weights)
            handle.write(
                "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (
                    holdout.family,
                    holdout.n_test,
                    holdout.n_vle,
                    holdout.k,
                    _fmt(holdout.train_point.recall),
                    _fmt(holdout.train_point.precision),
                    _fmt(holdout.train_point.threshold),
                    _fmt(_rate(holdout.holdout_called, holdout.n_vle)),
                    _fmt(_rate(holdout.stem_called, holdout.n_vle)),
                    weights,
                )
            )


def _summary_text(
    n_vle: int,
    n_decoy: int,
    terms: list[tuple[str, str]],
    published: _Model,
    search: _Search,
    weight_rows: list[tuple[str, str, int, float]],
    holdouts: list[_Outer],
) -> str:
    voting = [item for item in holdouts if item.voting]
    mean_recall = 0.0
    mean_stem = 0.0
    if voting:
        mean_recall = sum(_rate(item.holdout_called, item.n_vle) for item in voting) / len(voting)
        mean_stem = sum(_rate(item.stem_called, item.n_vle) for item in voting) / len(voting)
    lines = [
        "VLE discriminator score.",
        "Objective: maximize family-weighted recall at precision >= %.2f." % _TARGET_PRECISION,
        "Each VLE family with at least %d circles has the same total weight." % _MIN_FAMILY_N,
        "Those totals, plus weight 1 on smaller VLE families, sum to the VLE row count.",
        "VLE families below %d keep weight 1 per row." % _MIN_FAMILY_N,
        "Each DECOY family has the same total weight, and those totals sum to the DECOY row count.",
        "Precision counts a decoy by its row weight. The %.1f%% VLE prevalence therefore stays inside the cutoff."
        % (100.0 * n_vle / (n_vle + n_decoy)),
        "z(x) = (x - train_median) / train_IQR. Counts and lengths are log1p before that z-score.",
        "Each training fold orients a feature before selection. The sign is -1 when -z has the better family-weighted recall at precision >= %.2f. A tie keeps +1."
        % _TARGET_PRECISION,
        "A feature is eligible when its training IQR is positive.",
        "%s is excluded. An exact duplicate column is dropped." % _LENGTH_FEATURE,
        "After a term is accepted, a remaining candidate is dropped when its family-weighted absolute Pearson correlation with that term is at least %.2f."
        % _CORR_DROP,
        "Forward selection adds the mix that most raises training recall, and stops when recall does not rise.",
        "For each voting family, k is chosen from the other voting families.",
        "A longer prefix is kept when the unweighted mean of the held-out recall gains exceeds its standard error.",
        "The published k is the median of those choices. An even median rounds down.",
        "The published formula is that prefix of the all-data path, then a short Hooke-Jeeves polish.",
        "Validation holds out one VLE family at a time. Decoy rows stay in every training set.",
        "Held-out recall is the raw call rate of that family. Families below %d are scored and do not vote."
        % _MIN_FAMILY_N,
        "Each family ROC ranks that family against every decoy. AUCs are in auc_summary.txt.",
        "The stem baseline is the oriented log1p(%s) column." % _STEM_FEATURE,
        "Weights are nonnegative and sum to 1 on the oriented features. Call VLE when S >= threshold.",
        "Leave-one-family weights are nonnegative weights on that fold's oriented features.",
        "The threshold row in parameters.tsv stores the cutoff in the weight column.",
        "",
        "classes",
        "VLE\t%d" % n_vle,
        "DECOY\t%d" % n_decoy,
        "voting_families\t%d" % len(voting),
        "published_k\t%d" % published.k,
        "median_k\t%d" % published.requested_k,
        "",
        "family weights",
    ]
    for kind, name, count, total in weight_rows:
        denom = n_vle if kind == "VLE" else n_decoy
        lines.append(
            "%s\t%s\tn %d\ttotal_weight %s\tshare %s"
            % (kind, name, count, _fmt(total), _fmt(total / denom if denom else 0.0))
        )
    lines.append("")
    lines.append("features")
    drop_reason = {name: reason for name, reason, _rho in search.drops}
    chosen = {int(index) for index in np.flatnonzero(published.weights > 0.0)}
    for index, (name, _transform) in enumerate(terms):
        if name in drop_reason:
            status = drop_reason[name]
        elif index in chosen:
            status = "selected"
        else:
            status = "not selected"
        lines.append(
            "%s\t%s\trho %s\tsign %+d"
            % (name, status, _fmt(float(search.rho[index])), int(search.sign[index]))
        )
    lines.append("")
    lines.append("all-data path")
    for index, step in enumerate(search.steps, start=1):
        lines.append(
            "k %d\t%s\ttrain_weighted_recall %s"
            % (index, terms[step.feature][0], _fmt(step.train_point.recall))
        )
    lines.append("")
    lines.append("formula")
    pieces = []
    for index, (name, transform) in enumerate(terms):
        if published.weights[index] <= 0.0:
            continue
        source = "log1p(%s)" % name if transform == "log1p" else name
        lines.append(
            "z%d = (%s - %s) / %s"
            % (
                index + 1,
                source,
                _fmt(float(published.median[index])),
                _fmt(float(published.iqr[index])),
            )
        )
        term = "z%d" % (index + 1)
        if published.sign[index] < 0.0:
            term = "(-%s)" % term
        pieces.append("%s*%s" % (_fmt(float(published.weights[index])), term))
    lines.append("S = " + " + ".join(pieces))
    lines.append("threshold = %s" % _fmt(published.train_point.threshold))
    lines.append(
        "all-data train weighted recall %.4f at precision %.4f"
        % (published.train_point.recall, published.train_point.precision)
    )
    lines.append("polish_evals\t%d" % published.n_eval)
    lines.append("")
    lines.append("complexity")
    for holdout in holdouts:
        if not holdout.voting:
            continue
        lines.append(
            "%s\tchosen_k %d\tpath\t%s" % (holdout.family, holdout.k, _feature_names(terms, holdout.path_features))
        )
        gain_by_k = {row[0]: row for row in holdout.gain_rows}
        shown = max([1] + [row[0] for row in holdout.gain_rows])
        for k in range(1, shown + 1):
            recalls = holdout.inner_recalls[k - 1]
            parts = []
            for name, recall in zip(holdout.inner_names, recalls):
                parts.append("%s %.4f" % (name, recall))
            if k == 1:
                lines.append("k %d\t%s" % (k, "\t".join(parts)))
                continue
            _k, mean, se, accepted = gain_by_k[k]
            lines.append(
                "k %d\t%s\tmean_gain %s\tse %s\tkeep %s"
                % (k, "\t".join(parts), _fmt(mean), _fmt(se), "yes" if accepted else "no")
            )
    lines.append("")
    lines.append(
        "voting-family mean holdout recall %.4f, stem recall %.4f" % (mean_recall, mean_stem)
    )
    lines.append("")
    lines.append("holdouts")
    for holdout in holdouts:
        lines.append(
            "%s\tn %d\tk %d\tvote %s\tholdout_recall %.4f\tstem_recall %.4f\ttrain_weighted_recall %.4f"
            % (
                holdout.family,
                holdout.n_vle,
                holdout.k,
                "yes" if holdout.voting else "no",
                _rate(holdout.holdout_called, holdout.n_vle),
                _rate(holdout.stem_called, holdout.n_vle),
                holdout.train_point.recall,
            )
        )
    lines.append("")
    return "\n".join(lines)


def _save_fig(fig: plt.Figure, path: Path) -> None:
    temporary = path.with_suffix(".tmp.png")
    fig.savefig(temporary, dpi=150)
    plt.close(fig)
    temporary.replace(path)


def _write_weights(path: Path, terms: list[tuple[str, str]], holdouts: list[_Outer]) -> None:
    matrix = np.vstack([item.weights for item in holdouts])
    names = [name for name, _kind in terms]
    families = [item.family for item in holdouts]
    fig_w = max(12.0, 0.32 * len(names) + 3.0)
    fig_h = max(6.0, 0.32 * len(families) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    image = ax.imshow(matrix, vmin=0.0, vmax=1.0, cmap="viridis", aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(families)))
    ax.set_yticklabels(families, fontsize=8)
    ax.set_title("Leave-one-VLE-family-out weights")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Weight")
    fig.tight_layout()
    _save_fig(fig, path)


def _write_family_rates(path: Path, holdouts: list[_Outer]) -> None:
    names = [item.family for item in holdouts]
    rates = [_rate(item.holdout_called, item.n_vle) for item in holdouts]
    positions = np.arange(len(names))
    fig_height = max(6.0, 0.32 * len(names) + 1.4)
    fig, ax = plt.subplots(figsize=(10.5, fig_height))
    ax.barh(positions, rates, color=_COLOR_VLE, height=0.72)
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("Held-out recall")
    ax.set_title("Leave-one-VLE-family-out recall")
    fig.tight_layout()
    _save_fig(fig, path)


def _write_roc(path: Path, holdouts: list[_Outer]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.plot([0.0, 1.0], [0.0, 1.0], color="#bbbbbb", linewidth=1.0, zorder=0)
    for holdout in holdouts:
        ax.plot(
            holdout.fpr,
            holdout.tpr,
            linewidth=1.6,
            label="%s  %.3f" % (holdout.family, holdout.auc),
        )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Leave-one-VLE-family-out ROC")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(title="AUC", frameon=False, fontsize=8)
    fig.tight_layout()
    _save_fig(fig, path)


def _auc_summary_text(holdouts: list[_Outer]) -> str:
    voting = [item for item in holdouts if item.voting]
    lines = [
        "Leave-one-VLE-family-out ROC.",
        "Positives are the held-out VLE family. Negatives are every decoy row.",
        "The score is that family's formula, refit without the family.",
        "AUC is the unweighted Mann-Whitney statistic. Ties use the average rank.",
        "stem_auc is the oriented log1p(%s) column, with the same positives and negatives." % _STEM_FEATURE,
        "",
        "family\tn\tk\tauc\tstem_auc",
    ]
    for holdout in holdouts:
        lines.append(
            "%s\t%d\t%d\t%.4f\t%.4f"
            % (holdout.family, holdout.n_vle, holdout.k, holdout.auc, holdout.stem_auc)
        )
    if voting:
        mean_auc = sum(item.auc for item in voting) / len(voting)
        mean_stem = sum(item.stem_auc for item in voting) / len(voting)
        lines.append("")
        lines.append("voting-family mean auc %.4f, stem auc %.4f" % (mean_auc, mean_stem))
    lines.append("")
    return "\n".join(lines)


def _publish(search: _Search, requested_k: int, raw: np.ndarray, is_vle: np.ndarray, families: np.ndarray) -> _Model:
    used = min(requested_k, len(search.steps))
    if used < 1:
        _die("ERROR: the published path has no terms")
    start = search.steps[used - 1].weights.copy()
    support = np.flatnonzero(start > 0.0)
    if support.size == 0:
        _die("ERROR: the published prefix has no positive weight")
    terms = _WORKER["terms"]
    z = _apply_z(raw, search.median, search.iqr, search.active, terms) * search.sign
    row_weight = _row_weights(families, is_vle)
    weights, point, n_eval = _polish(z, is_vle, row_weight, start, support)
    print(
        "polished k %d  recall %.4f  evals %d" % (used, point.recall, n_eval),
        file=sys.stderr,
        flush=True,
    )
    return _Model(weights, search.median, search.iqr, point, n_eval, used, requested_k, search.sign)


def main() -> None:
    if len(sys.argv) > 3:
        _die("ERROR: usage: construct_VLE_discriminator_score.py [input.tsv] [output_dir]")
    _self_check()
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else _DEFAULT_INPUT
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else _DEFAULT_OUTPUT
    started = time.perf_counter()
    raw, is_vle, families, terms = _load(input_path)
    n_vle = int(is_vle.sum())
    n_decoy = int(is_vle.size - n_vle)
    try:
        stem_index = next(index for index, (name, _kind) in enumerate(terms) if name == _STEM_FEATURE)
    except StopIteration:
        _die("ERROR: missing %s" % _STEM_FEATURE)
    try:
        length_index = next(index for index, (name, _kind) in enumerate(terms) if name == _LENGTH_FEATURE)
    except StopIteration:
        _die("ERROR: missing %s" % _LENGTH_FEATURE)
    vle_names = _vle_family_names(families, is_vle)
    voting = _voting_names(families, is_vle)
    if len(voting) < 3:
        _die("ERROR: need at least 3 VLE families with %d or more circles" % _MIN_FAMILY_N)
    family_rows = {name: np.flatnonzero(families == name) for name in vle_names}
    print(
        "loaded VLE %d DECOY %d VLE-families %d voting %d in %.1fs"
        % (n_vle, n_decoy, len(vle_names), len(voting), time.perf_counter() - started),
        flush=True,
    )
    initargs = (raw, is_vle, families, terms, stem_index, length_index, voting, family_rows)
    _init_worker(*initargs)
    jobs = [("outer", family) for family in voting] + [("path", "")]
    workers = min(len(jobs), multiprocessing.cpu_count())
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(processes=workers, initializer=_init_worker, initargs=initargs) as pool:
        results = pool.map(_run_job, jobs)
    outers: list[_Outer] = []
    search: _Search | None = None
    for kind, payload in results:
        if kind == "path":
            search = payload
        else:
            outers.append(payload)
    if search is None or len(outers) != len(voting):
        _die("ERROR: a fit did not return")
    outers.sort(key=lambda item: item.family)
    for holdout in outers:
        print(
            "family %s  k %d  holdout_recall %.4f  auc %.4f  stem_recall %.4f  stem_auc %.4f"
            % (
                holdout.family,
                holdout.k,
                _rate(holdout.holdout_called, holdout.n_vle),
                holdout.auc,
                _rate(holdout.stem_called, holdout.n_vle),
                holdout.stem_auc,
            ),
            flush=True,
        )
    requested_k = _median_k([item.k for item in outers])
    published = _publish(search, requested_k, raw, is_vle, families)
    feature_order = tuple(step.feature for step in search.steps)
    holdouts = list(outers)
    for name in vle_names:
        if name in voting:
            continue
        holdouts.append(_holdout_fixed(name, feature_order, published.k))
    holdouts.sort(key=lambda item: item.family)
    for holdout in holdouts:
        if holdout.voting:
            continue
        print(
            "family %s  k %d  holdout_recall %.4f  auc %.4f  stem_recall %.4f  stem_auc %.4f"
            % (
                holdout.family,
                holdout.k,
                _rate(holdout.holdout_called, holdout.n_vle),
                holdout.auc,
                _rate(holdout.stem_called, holdout.n_vle),
                holdout.stem_auc,
            ),
            flush=True,
        )
    weight_rows = []
    all_weights = _row_weights(families, is_vle)
    for name in vle_names:
        member = (families == name) & is_vle
        weight_rows.append(("VLE", name, int(member.sum()), float(all_weights[member].sum())))
    decoy_names = sorted({str(name) for name in families[~is_vle]})
    for name in decoy_names:
        member = (families == name) & ~is_vle
        weight_rows.append(("DECOY", name, int(member.sum()), float(all_weights[member].sum())))
    output_dir.mkdir(parents=True, exist_ok=True)
    for stale in ("score_distributions.png", "precision_recall.png", "roc.png"):
        old = output_dir / stale
        if old.exists():
            old.unlink()
    summary_path = output_dir / "summary.txt"
    parameters_path = output_dir / "parameters.tsv"
    holdouts_path = output_dir / "holdouts.tsv"
    weights_path = output_dir / "weights.png"
    rates_path = output_dir / "family_rates.png"
    roc_path = output_dir / "holdout_roc.png"
    auc_path = output_dir / "auc_summary.txt"
    _write_parameters(parameters_path, terms, published)
    _write_holdouts(holdouts_path, terms, holdouts)
    text = _summary_text(n_vle, n_decoy, terms, published, search, weight_rows, holdouts)
    summary_path.write_text(text, encoding="utf-8", newline="\n")
    auc_text = _auc_summary_text(holdouts)
    auc_path.write_text(auc_text, encoding="utf-8", newline="\n")
    _write_weights(weights_path, terms, holdouts)
    _write_family_rates(rates_path, holdouts)
    _write_roc(roc_path, holdouts)
    print(text, end="")
    print(auc_text, end="")
    for path in (summary_path, parameters_path, holdouts_path, weights_path, rates_path, roc_path, auc_path):
        print("wrote %s" % path)
    print("elapsed %.1fs" % (time.perf_counter() - started), flush=True)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
