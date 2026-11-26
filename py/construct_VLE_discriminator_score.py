#!/usr/bin/env python3
"""Fit a closed-form VLE score and a precision-0.5 operating point.

S is a nonnegative weighted sum of four robust z-scores. Hooke-Jeeves
maximizes recall at precision >= 0.5 on a training half. The published
formula is split 0. Held-out figures are the precision-recall curve
(precision >= 0.25), the ROC with AUC, the score densities, and the weights.
"""

from __future__ import annotations

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
_N_SPLITS = 5
_SEED = 1
_TARGET_PRECISION = 0.5
_PR_FLOOR = 0.25
_MAX_EVALS = 400
_COLOR_VLE = "#d55e00"
_COLOR_DECOY = "#0072b2"
_COLOR_STEM = "#0072b2"

# w1..w4 in this order. The stem count is log1p before the z-score.
_TERMS = (
    ("longest_composite_stem_total_nr_bases", "log1p"),
    ("longest_composite_stem_pct_total_bases", "identity"),
    ("pct_bulge", "identity"),
    ("pct_paired_bases", "identity"),
)
_SHORT = (
    "log1p(composite stem)",
    "composite stem %",
    "pct_bulge",
    "pct_paired",
)


class _Point(NamedTuple):
    recall: float
    precision: float
    fpr: float
    threshold: float


class _SplitResult(NamedTuple):
    weights: np.ndarray
    median: np.ndarray
    iqr: np.ndarray
    train_point: _Point
    test_point: _Point
    deployed: _Point
    test_auc: float
    stem_test_point: _Point
    stem_test_auc: float
    n_eval: int
    test_scores: np.ndarray | None
    test_is_vle: np.ndarray | None
    stem_test_scores: np.ndarray | None


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


def _class_is_vle(label: str, line_number: int) -> bool:
    prefix = label.split("_", 1)[0]
    if prefix == "VLE":
        return True
    if prefix == "DECOY":
        return False
    _die("ERROR: line %d label '%s' is not VLE or DECOY" % (line_number, label))


def _load(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return the four raw features, shape (n, 4), and a VLE mask."""
    wanted = [name for name, _transform in _TERMS]
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
            raw_columns: list[list[float]] = [[], [], [], []]
            is_vle_list: list[bool] = []

            def consume(line_number: int, fields: list[str]) -> None:
                if len(fields) != width:
                    _die(
                        "ERROR: %s:%d has %d columns, expected %d"
                        % (path, line_number, len(fields), width)
                    )
                is_vle_list.append(_class_is_vle(fields[0], line_number))
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
    if not bool(is_vle.any()) or not bool((~is_vle).any()):
        _die("ERROR: need both VLE and DECOY rows")
    if np.any(raw[:, 0] < 0.0):
        _die("ERROR: longest_composite_stem_total_nr_bases is negative")
    return raw, is_vle


def _stratified_split(is_vle: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    train_parts: list[np.ndarray] = []
    test_parts: list[np.ndarray] = []
    for class_value in (True, False):
        members = np.flatnonzero(is_vle == class_value)
        members = members[rng.permutation(members.size)]
        n_test = members.size // 2
        test_parts.append(members[:n_test])
        train_parts.append(members[n_test:])
    return np.concatenate(train_parts), np.concatenate(test_parts)


def _transform(raw: np.ndarray) -> np.ndarray:
    transformed = np.empty_like(raw)
    transformed[:, 0] = np.log1p(raw[:, 0])
    transformed[:, 1:] = raw[:, 1:]
    return transformed


def _fit_z(train_raw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    transformed = _transform(train_raw)
    median = np.median(transformed, axis=0)
    q25, q75 = np.quantile(transformed, [0.25, 0.75], axis=0)
    iqr = q75 - q25
    if np.any(iqr <= 0.0):
        _die("ERROR: a training feature has IQR 0")
    return median, iqr


def _apply_z(raw: np.ndarray, median: np.ndarray, iqr: np.ndarray) -> np.ndarray:
    return (_transform(raw) - median) / iqr


def _tie_ends(sorted_scores: np.ndarray) -> np.ndarray:
    """Last index of each tied group in a descending score array."""
    change = np.empty(sorted_scores.size, dtype=bool)
    change[-1] = True
    change[:-1] = sorted_scores[:-1] != sorted_scores[1:]
    return np.flatnonzero(change)


def _ranked_curves(
    scores: np.ndarray, is_vle: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Tie-end recall, precision, FPR, TPR, and the score of each group."""
    n_pos = int(is_vle.sum())
    n_neg = int(scores.size - n_pos)
    if n_pos == 0 or n_neg == 0:
        _die("ERROR: a split is missing a class")
    order = np.argsort(-scores, kind="mergesort")
    sorted_scores = scores[order]
    true_positive = np.cumsum(is_vle[order], dtype=np.int64)
    ends = _tie_ends(sorted_scores)
    tp = true_positive[ends].astype(np.float64)
    predicted = (ends + 1).astype(np.float64)
    fp = predicted - tp
    recall = tp / n_pos
    precision = tp / predicted
    fpr = fp / n_neg
    return recall, precision, fpr, recall.copy(), sorted_scores[ends]


def _point_at_precision(
    scores: np.ndarray, is_vle: np.ndarray, target: float
) -> _Point:
    """Largest recall whose threshold still has precision >= target.

    Ties share one threshold. Among equal recalls, the higher threshold wins.
    """
    recall, precision, fpr, _tpr, group_score = _ranked_curves(scores, is_vle)
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


def _deployed_point(scores: np.ndarray, is_vle: np.ndarray, threshold: float) -> _Point:
    if not np.isfinite(threshold):
        return _Point(0.0, 0.0, 0.0, threshold)
    n_pos = int(is_vle.sum())
    n_neg = int(scores.size - n_pos)
    called = scores >= threshold
    tp = float(np.count_nonzero(called & is_vle))
    fp = float(np.count_nonzero(called & ~is_vle))
    predicted = tp + fp
    precision = tp / predicted if predicted else 0.0
    recall = tp / n_pos if n_pos else 0.0
    fpr = fp / n_neg if n_neg else 0.0
    return _Point(recall, precision, fpr, threshold)


def _auc(scores: np.ndarray, is_vle: np.ndarray) -> float:
    n_pos = int(is_vle.sum())
    n_neg = int(scores.size - n_pos)
    if n_pos == 0 or n_neg == 0 or np.all(scores == scores[0]):
        return 0.5
    order = np.argsort(scores, kind="mergesort")
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


def _project(weights: np.ndarray) -> np.ndarray:
    clipped = np.maximum(weights, 0.0)
    total = float(clipped.sum())
    if total <= 0.0:
        return np.full(weights.size, 1.0 / weights.size)
    return clipped / total


def _fisher_start(features: np.ndarray, is_vle: np.ndarray) -> np.ndarray:
    """Equal-class Fisher direction, projected onto the simplex."""
    positive = features[is_vle]
    negative = features[~is_vle]
    mean_gap = positive.mean(axis=0) - negative.mean(axis=0)

    def _covariance(block: np.ndarray) -> np.ndarray:
        centered = block - block.mean(axis=0)
        return (centered.T @ centered) / max(block.shape[0] - 1, 1)

    sigma = 0.5 * (_covariance(positive) + _covariance(negative))
    direction, _, _, _ = np.linalg.lstsq(sigma, mean_gap, rcond=None)
    if float(np.dot(direction, mean_gap)) < 0.0:
        direction = -direction
    return _project(direction)


class _Budget:
    def __init__(self, limit: int):
        self.limit = limit
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
        point = _point_at_precision(features @ weights, is_vle, _TARGET_PRECISION)
        self.cache[key] = point.recall
        return point.recall


def _exploratory(
    origin: np.ndarray,
    origin_value: float,
    delta: float,
    budget: _Budget,
    features: np.ndarray,
    is_vle: np.ndarray,
) -> tuple[np.ndarray, float]:
    current = origin.copy()
    current_value = origin_value
    for coordinate in range(current.size):
        best = current
        best_value = current_value
        for sign in (1.0, -1.0):
            candidate = current.copy()
            candidate[coordinate] += sign * delta
            candidate = _project(candidate)
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
) -> tuple[np.ndarray, float]:
    base = _project(start)
    base_value = budget.recall(base, features, is_vle)
    if base_value is None:
        _die("ERROR: Hooke-Jeeves has no evaluations left at the start")
    delta = 0.25
    while delta >= 0.01 and not budget.exhausted:
        trial, trial_value = _exploratory(base, base_value, delta, budget, features, is_vle)
        if trial_value <= base_value + 1e-12:
            delta *= 0.5
            continue
        previous = base
        base = trial
        base_value = trial_value
        jumps = 0
        while jumps < 25 and not budget.exhausted:
            jumps += 1
            pattern = _project(base + (base - previous))
            pattern_value = budget.recall(pattern, features, is_vle)
            if pattern_value is None:
                break
            anchor = pattern_value if pattern_value > base_value else base_value
            refined, refined_value = _exploratory(
                pattern, anchor, delta, budget, features, is_vle
            )
            if refined_value > base_value + 1e-12:
                previous = base
                base = refined
                base_value = refined_value
            else:
                break
    return base, base_value


def _fit_weights(features: np.ndarray, is_vle: np.ndarray) -> tuple[np.ndarray, float, int]:
    cache: dict[tuple[float, ...], float] = {}
    starts = [np.eye(features.shape[1])[index] for index in range(features.shape[1])]
    starts.append(np.full(features.shape[1], 1.0 / features.shape[1]))
    starts.append(_fisher_start(features, is_vle))
    best_weights = _project(starts[0])
    best_recall = -1.0
    best_auc = -1.0
    total_eval = 0
    for index, start in enumerate(starts):
        budget = _Budget(_MAX_EVALS)
        budget.cache = cache
        weights, recall = _hooke_jeeves(start, budget, features, is_vle)
        total_eval += budget.n_eval
        auc = _auc(features @ weights, is_vle)
        if recall > best_recall + 1e-12 or (abs(recall - best_recall) <= 1e-12 and auc > best_auc):
            best_weights = weights
            best_recall = recall
            best_auc = auc
        print(
            "  start %d  recall@P>=%.2f %.4f  auc %.4f  new_evals %d"
            % (index + 1, _TARGET_PRECISION, recall, auc, budget.n_eval),
            file=sys.stderr,
            flush=True,
        )
    return best_weights, best_recall, total_eval


def _fit_split(
    raw: np.ndarray,
    is_vle: np.ndarray,
    train_index: np.ndarray,
    test_index: np.ndarray,
    keep_scores: bool,
) -> _SplitResult:
    median, iqr = _fit_z(raw[train_index])
    train_z = _apply_z(raw[train_index], median, iqr)
    test_z = _apply_z(raw[test_index], median, iqr)
    train_class = is_vle[train_index]
    test_class = is_vle[test_index]
    weights, _train_recall, n_eval = _fit_weights(train_z, train_class)
    train_scores = train_z @ weights
    test_scores = test_z @ weights
    stem_test = test_z[:, 0]
    train_point = _point_at_precision(train_scores, train_class, _TARGET_PRECISION)
    test_point = _point_at_precision(test_scores, test_class, _TARGET_PRECISION)
    deployed = _deployed_point(test_scores, test_class, train_point.threshold)
    stem_test_point = _point_at_precision(stem_test, test_class, _TARGET_PRECISION)
    return _SplitResult(
        weights=weights,
        median=median,
        iqr=iqr,
        train_point=train_point,
        test_point=test_point,
        deployed=deployed,
        test_auc=_auc(test_scores, test_class),
        stem_test_point=stem_test_point,
        stem_test_auc=_auc(stem_test, test_class),
        n_eval=n_eval,
        test_scores=test_scores if keep_scores else None,
        test_is_vle=test_class if keep_scores else None,
        stem_test_scores=stem_test if keep_scores else None,
    )


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
    projected = _project(np.array([1.0, -2.0, 1.0, 0.0]))
    _expect(abs(float(projected.sum()) - 1.0) < 1e-12, "simplex sums to 1")
    _expect(float(projected.min()) >= -1e-12, "simplex is nonnegative")


def _write_parameters(path: Path, result: _SplitResult) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("feature\ttransform\ttrain_median\ttrain_iqr\tweight\n")
        for index, (name, transform) in enumerate(_TERMS):
            handle.write(
                "%s\t%s\t%s\t%s\t%s\n"
                % (
                    name,
                    transform,
                    _fmt(float(result.median[index])),
                    _fmt(float(result.iqr[index])),
                    _fmt(float(result.weights[index])),
                )
            )
        handle.write("threshold\t\t\t\t%s\n" % _fmt(result.train_point.threshold))


def _write_holdouts(path: Path, results: list[_SplitResult]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "split\tw1\tw2\tw3\tw4\t"
            "train_recall_at_p50\ttrain_precision\t"
            "test_recall_at_p50\ttest_precision_at_p50\ttest_auc\t"
            "stem_test_recall_at_p50\tstem_test_precision_at_p50\tstem_test_auc\t"
            "deployed_precision\tdeployed_recall\tdeployed_fpr\ttrain_threshold\n"
        )
        for split_id, result in enumerate(results):
            weights = "\t".join(_fmt(float(value)) for value in result.weights)
            handle.write(
                "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (
                    split_id,
                    weights,
                    _fmt(result.train_point.recall),
                    _fmt(result.train_point.precision),
                    _fmt(result.test_point.recall),
                    _fmt(result.test_point.precision),
                    _fmt(result.test_auc),
                    _fmt(result.stem_test_point.recall),
                    _fmt(result.stem_test_point.precision),
                    _fmt(result.stem_test_auc),
                    _fmt(result.deployed.precision),
                    _fmt(result.deployed.recall),
                    _fmt(result.deployed.fpr),
                    _fmt(result.train_point.threshold),
                )
            )


def _summary_text(n_vle: int, n_decoy: int, results: list[_SplitResult]) -> str:
    published = results[0]
    lines = [
        "VLE discriminator score.",
        "Objective: maximize recall at precision >= %.2f." % _TARGET_PRECISION,
        "Precision and recall use raw row counts, so the %.1f%% VLE prevalence is part of the objective."
        % (100.0 * n_vle / (n_vle + n_decoy)),
        "z(x) = (x - train_median) / train_IQR. The stem count is log1p before that z-score.",
        "Weights are nonnegative and sum to 1. Call VLE when S >= threshold.",
        "The threshold row in parameters.tsv stores the cutoff in the weight column.",
        "test_recall_at_p50 chooses the threshold on the held-out half.",
        "deployed_* applies the training threshold to the held-out half.",
        "The stem columns are w = (1, 0, 0, 0): log1p(composite stem) alone.",
        "",
        "classes",
        "VLE\t%d" % n_vle,
        "DECOY\t%d" % n_decoy,
        "",
        "published split\t0",
        "formula",
    ]
    terms = []
    for index, (name, transform) in enumerate(_TERMS):
        source = "log1p(%s)" % name if transform == "log1p" else name
        lines.append(
            "z%d = (%s - %s) / %s"
            % (index + 1, source, _fmt(float(published.median[index])), _fmt(float(published.iqr[index])))
        )
        terms.append("%s*z%d" % (_fmt(float(published.weights[index])), index + 1))
    lines.append("S = " + " + ".join(terms))
    lines.append("threshold = %s" % _fmt(published.train_point.threshold))
    lines.append(
        "split 0 train recall %.4f at precision %.4f"
        % (published.train_point.recall, published.train_point.precision)
    )
    lines.append(
        "split 0 test recall %.4f at precision %.4f, AUC %.4f"
        % (published.test_point.recall, published.test_point.precision, published.test_auc)
    )
    lines.append(
        "split 0 deployed precision %.4f, recall %.4f, FPR %.6f"
        % (published.deployed.precision, published.deployed.recall, published.deployed.fpr)
    )
    lines.append(
        "split 0 log-stem alone test recall %.4f at precision %.4f, AUC %.4f"
        % (
            published.stem_test_point.recall,
            published.stem_test_point.precision,
            published.stem_test_auc,
        )
    )
    lines.append("")
    lines.append("holdouts")
    for split_id, result in enumerate(results):
        weight_text = ", ".join(_fmt(float(value)) for value in result.weights)
        lines.append(
            "%d\tweights %s\ttrain_recall %.4f\ttest_recall %.4f\ttest_auc %.4f\tstem_test_recall %.4f\tstem_test_auc %.4f"
            % (
                split_id,
                weight_text,
                result.train_point.recall,
                result.test_point.recall,
                result.test_auc,
                result.stem_test_point.recall,
                result.stem_test_auc,
            )
        )
    lines.append("")
    return "\n".join(lines)


def _save_fig(fig: plt.Figure, path: Path) -> None:
    temporary = path.with_suffix(".tmp.png")
    fig.savefig(temporary, dpi=150)
    plt.close(fig)
    temporary.replace(path)


def _write_distributions(path: Path, scores: np.ndarray, is_vle: np.ndarray, threshold: float) -> None:
    low, high = np.quantile(scores, [0.005, 0.995])
    low = float(low)
    high = float(high)
    if np.isfinite(threshold):
        low = min(low, threshold)
        high = max(high, threshold)
    if not high > low:
        low -= 0.5
        high += 0.5
    bins = np.linspace(low, high, 81)
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    decoy = scores[~is_vle]
    vle = scores[is_vle]
    decoy = decoy[(decoy >= low) & (decoy <= high)]
    vle = vle[(vle >= low) & (vle <= high)]
    if decoy.size:
        ax.hist(decoy, bins=bins, density=True, color=_COLOR_DECOY, alpha=0.55, label="DECOY")
    if vle.size:
        ax.hist(vle, bins=bins, density=True, color=_COLOR_VLE, alpha=0.55, label="VLE")
    if np.isfinite(threshold):
        ax.axvline(threshold, color="0.1", linewidth=1.2, label="train threshold")
    ax.set_xlim(low, high)
    ax.set_xlabel("Score S")
    ax.set_ylabel("Density")
    ax.set_title("Held-out score distribution, split 0")
    ax.legend(frameon=False)
    fig.tight_layout()
    _save_fig(fig, path)


def _write_pr(
    path: Path,
    scores: np.ndarray,
    stem_scores: np.ndarray,
    is_vle: np.ndarray,
    deployed: _Point,
    score_point: _Point,
    stem_point: _Point,
) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    for series, color, label, point in (
        (scores, _COLOR_VLE, "S", score_point),
        (stem_scores, _COLOR_STEM, "log stem", stem_point),
    ):
        recall, precision, _fpr, _tpr, _group_score = _ranked_curves(series, is_vle)
        shown = precision.astype(np.float64).copy()
        shown[shown < _PR_FLOOR] = np.nan
        ax.plot(
            recall,
            shown,
            color=color,
            linewidth=1.8,
            label="%s  recall@P>=%.2f %.3f" % (label, _TARGET_PRECISION, point.recall),
        )
    if np.isfinite(deployed.threshold):
        ax.scatter(
            [deployed.recall],
            [deployed.precision],
            s=36,
            color="0.1",
            zorder=3,
            label="train threshold on test",
        )
    ax.axhline(_TARGET_PRECISION, color="0.45", linestyle="--", linewidth=1, label="precision %.2f" % _TARGET_PRECISION)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(_PR_FLOOR, 1.02)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Held-out precision-recall, precision >= %.2f" % _PR_FLOOR)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    _save_fig(fig, path)


def _write_roc(
    path: Path,
    scores: np.ndarray,
    stem_scores: np.ndarray,
    is_vle: np.ndarray,
    deployed: _Point,
    score_auc: float,
    stem_auc: float,
) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    ax.plot([0.0, 1.0], [0.0, 1.0], color="0.5", linestyle="--", linewidth=1, label="chance")
    for series, color, label, auc in (
        (scores, _COLOR_VLE, "S", score_auc),
        (stem_scores, _COLOR_STEM, "log stem", stem_auc),
    ):
        _recall, _precision, fpr, tpr, _group_score = _ranked_curves(series, is_vle)
        fpr = np.concatenate([[0.0], fpr])
        tpr = np.concatenate([[0.0], tpr])
        ax.plot(fpr, tpr, color=color, linewidth=1.8, label="%s  AUC %.3f" % (label, auc))
    if np.isfinite(deployed.threshold):
        ax.scatter(
            [deployed.fpr],
            [deployed.recall],
            s=36,
            color="0.1",
            zorder=3,
            label="train threshold on test",
        )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.02)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("Recall")
    ax.set_title("Held-out ROC")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    _save_fig(fig, path)


def _write_weights(path: Path, results: list[_SplitResult]) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    x = np.arange(len(_TERMS))
    width = 0.15
    for split_id, result in enumerate(results):
        ax.bar(
            x + (split_id - (_N_SPLITS - 1) / 2.0) * width,
            result.weights,
            width=width,
            label="split %d" % split_id,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(_SHORT, rotation=15, ha="right")
    ax.set_ylabel("Weight")
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Hooke-Jeeves weights by split")
    ax.legend(frameon=False, ncol=5, loc="upper right")
    fig.tight_layout()
    _save_fig(fig, path)


def main() -> None:
    if len(sys.argv) > 3:
        _die("ERROR: usage: construct_VLE_discriminator_score.py [input.tsv] [output_dir]")
    _self_check()
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else _DEFAULT_INPUT
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else _DEFAULT_OUTPUT
    started = time.perf_counter()
    raw, is_vle = _load(input_path)
    n_vle = int(is_vle.sum())
    n_decoy = int(is_vle.size - n_vle)
    print("loaded VLE %d DECOY %d in %.1fs" % (n_vle, n_decoy, time.perf_counter() - started), flush=True)
    rng = np.random.default_rng(_SEED)
    results: list[_SplitResult] = []
    for split_id in range(_N_SPLITS):
        train_index, test_index = _stratified_split(is_vle, rng)
        print("split %d  train %d  test %d" % (split_id, train_index.size, test_index.size), file=sys.stderr, flush=True)
        result = _fit_split(raw, is_vle, train_index, test_index, keep_scores=(split_id == 0))
        results.append(result)
        print(
            "split %d  test recall@P>=%.2f %.4f  AUC %.4f  stem recall %.4f  stem AUC %.4f"
            % (
                split_id,
                _TARGET_PRECISION,
                result.test_point.recall,
                result.test_auc,
                result.stem_test_point.recall,
                result.stem_test_auc,
            ),
            flush=True,
        )
    published = results[0]
    if published.test_scores is None or published.test_is_vle is None or published.stem_test_scores is None:
        _die("ERROR: split 0 did not keep held-out scores")
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "summary.txt"
    parameters_path = output_dir / "parameters.tsv"
    holdouts_path = output_dir / "holdouts.tsv"
    _write_parameters(parameters_path, published)
    _write_holdouts(holdouts_path, results)
    text = _summary_text(n_vle, n_decoy, results)
    summary_path.write_text(text, encoding="utf-8", newline="\n")
    _write_distributions(
        output_dir / "score_distributions.png",
        published.test_scores,
        published.test_is_vle,
        published.train_point.threshold,
    )
    _write_pr(
        output_dir / "precision_recall.png",
        published.test_scores,
        published.stem_test_scores,
        published.test_is_vle,
        published.deployed,
        published.test_point,
        published.stem_test_point,
    )
    _write_roc(
        output_dir / "roc.png",
        published.test_scores,
        published.stem_test_scores,
        published.test_is_vle,
        published.deployed,
        published.test_auc,
        published.stem_test_auc,
    )
    _write_weights(output_dir / "weights.png", results)
    print(text, end="")
    for path in (
        summary_path,
        parameters_path,
        holdouts_path,
        output_dir / "score_distributions.png",
        output_dir / "precision_recall.png",
        output_dir / "roc.png",
        output_dir / "weights.png",
    ):
        print("wrote %s" % path)
    print("elapsed %.1fs" % (time.perf_counter() - started), flush=True)


if __name__ == "__main__":
    main()
