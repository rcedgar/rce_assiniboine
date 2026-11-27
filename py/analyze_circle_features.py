#!/usr/bin/env python3
"""Rank circle-structure features by how well they separate VLE from DECOY.

Each row is one circle. The class is the label prefix (VLE or DECOY).
The columns are the combined ROD v2 structure features and RCE circle
features returned by extract_circular_rna_features. They are scored by a
tie-aware Mann-Whitney AUC and by Cohen's d. A second AUC is computed
after residualizing each feature on log(nr_bases), because raw counts
largely track circle length.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import NamedTuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_INPUT = (
    _ROOT / "big_output" / "controls_folds_100_per_family_circle_features_VLE_excludes_discarded.tsv"
)
_DEFAULT_OUTPUT = _ROOT / "analyze_circle_features"
_LENGTH_FEATURE = "nr_bases"
_NEAR_CONSTANT_FRACTION = 0.99
_TOP_N = 8
_COLOR_VLE = "#d55e00"
_COLOR_DECOY = "#0072b2"
_COLOR_TIE = "0.55"


class _Score(NamedTuple):
    name: str
    auc: float
    oriented: float
    higher: str
    cohens_d: float
    auc_adj: float | None
    oriented_adj: float | None
    higher_adj: str | None
    median_vle: float
    median_decoy: float
    mean_vle: float
    mean_decoy: float
    std_vle: float
    std_decoy: float


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


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
    _die(f"ERROR: line {line_number} label '{label}' is not VLE or DECOY")


def _family_of(label: str, line_number: int) -> str:
    parts = label.split("_", 2)
    if len(parts) != 3 or parts[0] not in ("VLE", "DECOY") or not parts[1] or not parts[2]:
        _die(
            "ERROR: line %d label '%s' is not VLE_family_... or DECOY_family_..."
            % (line_number, label)
        )
    return parts[1]


def _load(path: Path) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """Return feature names, values shaped (n, p), a VLE mask, and family names."""
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            first = handle.readline()
            if first == "":
                _die(f"ERROR: {path} is empty")
            first_fields = first.lstrip("\ufeff").rstrip("\r\n").split("\t")
            pending: list[tuple[int, list[str]]] = []
            if _is_header(first_fields[0]):
                feature_names = first_fields[1:]
            else:
                feature_names = _canonical_feature_names()
                pending.append((1, first_fields))
            if not feature_names:
                _die(f"ERROR: {path} has no feature columns")
            if len(feature_names) != len(set(feature_names)):
                _die(f"ERROR: {path} repeats a feature name")
            if _LENGTH_FEATURE not in feature_names:
                _die(f"ERROR: {path} has no {_LENGTH_FEATURE} column")
            width = 1 + len(feature_names)
            for line_number, line in enumerate(handle, start=2):
                text = line.rstrip("\r\n")
                if text == "":
                    continue
                pending.append((line_number, text.split("\t")))
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")

    if not pending:
        _die(f"ERROR: no data rows in {path}")

    values = np.empty((len(pending), len(feature_names)), dtype=np.float64)
    is_vle = np.empty(len(pending), dtype=bool)
    families = np.empty(len(pending), dtype=object)
    for row_index, (line_number, fields) in enumerate(pending):
        if len(fields) != width:
            _die(
                f"ERROR: {path}:{line_number} has {len(fields)} columns, expected {width}"
            )
        is_vle[row_index] = _class_is_vle(fields[0], line_number)
        families[row_index] = _family_of(fields[0], line_number)
        for feature_index, raw in enumerate(fields[1:]):
            try:
                value = float(raw)
            except ValueError:
                _die(
                    f"ERROR: {path}:{line_number} {feature_names[feature_index]} is not a number"
                )
            if not math.isfinite(value):
                _die(
                    f"ERROR: {path}:{line_number} {feature_names[feature_index]} is not finite"
                )
            values[row_index, feature_index] = value
    if not bool(is_vle.any()) or not bool((~is_vle).any()):
        _die("ERROR: need both VLE and DECOY rows")
    return feature_names, values, is_vle, families.astype(str)


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


def _auc(values: np.ndarray, is_vle: np.ndarray) -> tuple[float, str]:
    """Mann-Whitney AUC for VLE versus DECOY, with average ranks for ties.

    The returned AUC is P(VLE > DECOY) + 0.5 P(tie). Above 0.5, VLE tends
    to have the larger value.
    """
    n_vle = int(is_vle.sum())
    n_decoy = int(values.size - n_vle)
    if np.all(values == values[0]):
        return 0.5, "tie"
    rank_sum = float(_average_ranks(values)[is_vle].sum())
    auc = (rank_sum - n_vle * (n_vle + 1) / 2.0) / (n_vle * n_decoy)
    if auc > 0.5:
        return auc, "VLE"
    if auc < 0.5:
        return auc, "DECOY"
    return auc, "tie"


def _cohens_d(vle: np.ndarray, decoy: np.ndarray) -> float:
    """(mean VLE - mean DECOY) / pooled SD. A constant column returns 0."""
    n1 = vle.size
    n2 = decoy.size
    mean_diff = float(vle.mean() - decoy.mean())
    var1 = float(vle.var(ddof=1)) if n1 > 1 else 0.0
    var2 = float(decoy.var(ddof=1)) if n2 > 1 else 0.0
    pooled = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
    if pooled <= 0.0:
        if mean_diff == 0.0:
            return 0.0
        return math.copysign(math.inf, mean_diff)
    return mean_diff / math.sqrt(pooled)


def _residualize(values: np.ndarray, log_length: np.ndarray) -> np.ndarray:
    design = np.column_stack((np.ones(values.size), log_length))
    coefficient, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coefficient


def _class_std(values: np.ndarray) -> float:
    if values.size < 2:
        return 0.0
    return float(values.std(ddof=1))


def _score_features(
    feature_names: list[str],
    values: np.ndarray,
    is_vle: np.ndarray,
) -> list[_Score]:
    length_index = feature_names.index(_LENGTH_FEATURE)
    length = values[:, length_index]
    if np.any(length <= 0.0):
        _die(f"ERROR: {_LENGTH_FEATURE} must be positive")
    log_length = np.log(length)
    scores: list[_Score] = []
    for feature_index, name in enumerate(feature_names):
        column = values[:, feature_index]
        vle = column[is_vle]
        decoy = column[~is_vle]
        auc, higher = _auc(column, is_vle)
        if name == _LENGTH_FEATURE:
            auc_adj = None
            oriented_adj = None
            higher_adj = None
        else:
            auc_adj, higher_adj = _auc(_residualize(column, log_length), is_vle)
            oriented_adj = max(auc_adj, 1.0 - auc_adj)
        scores.append(
            _Score(
                name=name,
                auc=auc,
                oriented=max(auc, 1.0 - auc),
                higher=higher,
                cohens_d=_cohens_d(vle, decoy),
                auc_adj=auc_adj,
                oriented_adj=oriented_adj,
                higher_adj=higher_adj,
                median_vle=float(np.median(vle)),
                median_decoy=float(np.median(decoy)),
                mean_vle=float(vle.mean()),
                mean_decoy=float(decoy.mean()),
                std_vle=_class_std(vle),
                std_decoy=_class_std(decoy),
            )
        )
    scores.sort(key=lambda score: (-score.oriented, -abs(score.cohens_d), score.name))
    return scores


def _fmt(value: float | None) -> str:
    if value is None:
        return "NA"
    if math.isinf(value):
        return "inf" if value > 0.0 else "-inf"
    return "%.6g" % value


def _write_ranking(path: Path, scores: list[_Score]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "feature\thigher_class\tauc\tauc_oriented\tcohens_d\t"
            "auc_length_adjusted\tauc_oriented_length_adjusted\t"
            "higher_class_length_adjusted\t"
            "median_vle\tmedian_decoy\tmean_vle\tmean_decoy\tstd_vle\tstd_decoy\n"
        )
        for score in scores:
            higher_adj = "NA" if score.higher_adj is None else score.higher_adj
            handle.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (
                    score.name,
                    score.higher,
                    _fmt(score.auc),
                    _fmt(score.oriented),
                    _fmt(score.cohens_d),
                    _fmt(score.auc_adj),
                    _fmt(score.oriented_adj),
                    higher_adj,
                    _fmt(score.median_vle),
                    _fmt(score.median_decoy),
                    _fmt(score.mean_vle),
                    _fmt(score.mean_decoy),
                    _fmt(score.std_vle),
                    _fmt(score.std_decoy),
                )
            )


def _modal_fraction(values: np.ndarray) -> tuple[float, float]:
    unique, counts = np.unique(values, return_counts=True)
    mode_index = int(np.argmax(counts))
    return float(unique[mode_index]), float(counts[mode_index] / values.size)


def _near_constant(
    feature_names: list[str], values: np.ndarray
) -> list[tuple[str, float, float]]:
    found: list[tuple[str, float, float]] = []
    for feature_index, name in enumerate(feature_names):
        mode, fraction = _modal_fraction(values[:, feature_index])
        if fraction >= _NEAR_CONSTANT_FRACTION:
            found.append((name, mode, fraction))
    found.sort(key=lambda item: (-item[2], item[0]))
    return found


def _ranking_lines(scores: list[_Score], adjusted: bool) -> list[str]:
    lines: list[str] = []
    shown = scores[:_TOP_N]
    for rank, score in enumerate(shown, start=1):
        if adjusted:
            if score.oriented_adj is None or score.higher_adj is None:
                continue
            lines.append(
                "%d\t%s\toriented_auc=%.4f\thigher=%s\tmedian_vle=%s\tmedian_decoy=%s"
                % (
                    rank,
                    score.name,
                    score.oriented_adj,
                    score.higher_adj,
                    _fmt(score.median_vle),
                    _fmt(score.median_decoy),
                )
            )
        else:
            lines.append(
                "%d\t%s\toriented_auc=%.4f\thigher=%s\tcohens_d=%s\tmedian_vle=%s\tmedian_decoy=%s"
                % (
                    rank,
                    score.name,
                    score.oriented,
                    score.higher,
                    _fmt(score.cohens_d),
                    _fmt(score.median_vle),
                    _fmt(score.median_decoy),
                )
            )
    return lines


def _summary_text(
    n_vle: int,
    n_decoy: int,
    scores: list[_Score],
    adjusted: list[_Score],
    near_constant: list[tuple[str, float, float]],
) -> str:
    lines = [
        "Circle features, VLE versus DECOY.",
        "Oriented AUC is max(AUC, 1-AUC). Chance is 0.5.",
        "higher is the class whose values tend to be larger.",
        "Cohen's d is (mean VLE - mean DECOY) / pooled SD.",
        "",
        "classes",
        "VLE\t%d" % n_vle,
        "DECOY\t%d" % n_decoy,
        "",
        "top %d by raw oriented AUC" % _TOP_N,
    ]
    lines.extend(_ranking_lines(scores, adjusted=False))
    lines.append("")
    lines.append(
        "top %d by oriented AUC after residualizing on log(%s)"
        % (_TOP_N, _LENGTH_FEATURE)
    )
    lines.append("%s is the covariate and is omitted." % _LENGTH_FEATURE)
    lines.extend(_ranking_lines(adjusted, adjusted=True))
    lines.append("")
    lines.append(
        "near-constant features (at least %.0f%% of rows share one value)"
        % (100.0 * _NEAR_CONSTANT_FRACTION)
    )
    if not near_constant:
        lines.append("none")
    else:
        for name, mode, fraction in near_constant:
            lines.append("%s\tvalue=%s\tfraction=%.4f" % (name, _fmt(mode), fraction))
    lines.append("")
    lines.append("Full table: feature_ranking.tsv")
    lines.append("")
    return "\n".join(lines)


def _write_summary(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8", newline="\n")


def _bar_color(higher: str) -> str:
    if higher == "VLE":
        return _COLOR_VLE
    if higher == "DECOY":
        return _COLOR_DECOY
    return _COLOR_TIE


def _save_fig(fig: plt.Figure, path: Path) -> None:
    temporary = path.with_suffix(".tmp.png")
    fig.savefig(temporary, dpi=150)
    plt.close(fig)
    temporary.replace(path)


def _bar_rows(
    scores: list[_Score], adjusted: bool
) -> tuple[list[str], list[float], list[str]]:
    """Names, oriented AUC, and higher class, worst feature first."""
    if adjusted:
        chosen = [
            score
            for score in scores
            if score.oriented_adj is not None and score.higher_adj is not None
        ]
        chosen.sort(
            key=lambda score: (-float(score.oriented_adj or 0.0), -abs(score.cohens_d), score.name)
        )
        names = [score.name for score in chosen]
        oriented = [float(score.oriented_adj or 0.0) for score in chosen]
        higher = [str(score.higher_adj) for score in chosen]
    else:
        names = [score.name for score in scores]
        oriented = [score.oriented for score in scores]
        higher = [score.higher for score in scores]
    names.reverse()
    oriented.reverse()
    higher.reverse()
    return names, oriented, higher


def _write_ranking_bars(
    path: Path,
    scores: list[_Score],
    title: str,
    adjusted: bool,
) -> None:
    names, oriented, higher = _bar_rows(scores, adjusted)
    if not names:
        _die(f"ERROR: no features to plot in {path.name}")
    positions = np.arange(len(names))
    fig_height = max(6.0, 0.32 * len(names) + 1.4)
    fig, ax = plt.subplots(figsize=(10.5, fig_height))
    colors = [_bar_color(label) for label in higher]
    ax.barh(positions, np.asarray(oriented) - 0.5, left=0.5, color=colors, height=0.72, zorder=2)
    ax.axvline(0.5, color="0.15", linewidth=1, zorder=3)
    for position, score, label in zip(positions, oriented, higher):
        ax.text(score + 0.008, position, "%s  %.3f" % (label, score), va="center", ha="left", fontsize=8)
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlim(0.45, 1.18)
    ax.set_xlabel("Oriented AUC")
    ax.set_title(title)
    ax.grid(axis="x", color="0.85", zorder=0)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    _save_fig(fig, path)


def _panel_limits(column: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(column, [0.01, 0.99])
    low = float(low)
    high = float(high)
    if high > low:
        return low, high
    pad = 0.5 if low == 0.0 else max(0.5, abs(low) * 0.05)
    return low - pad, high + pad


def _write_distributions(
    path: Path,
    scores: list[_Score],
    feature_names: list[str],
    values: np.ndarray,
    is_vle: np.ndarray,
) -> None:
    n_features = len(scores)
    ncols = 6 if n_features > 6 else max(1, n_features)
    nrows = math.ceil(n_features / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.15 * ncols, 2.45 * nrows),
        squeeze=False,
    )
    flat = axes.ravel()
    name_to_index = {name: index for index, name in enumerate(feature_names)}
    for panel, score in enumerate(scores):
        ax = flat[panel]
        column = values[:, name_to_index[score.name]]
        low, high = _panel_limits(column)
        bins = np.linspace(low, high, 41)
        decoy = column[~is_vle]
        vle = column[is_vle]
        decoy = decoy[(decoy >= low) & (decoy <= high)]
        vle = vle[(vle >= low) & (vle <= high)]
        if decoy.size:
            ax.hist(decoy, bins=bins, density=True, color=_COLOR_DECOY, alpha=0.55, label="DECOY")
        if vle.size:
            ax.hist(vle, bins=bins, density=True, color=_COLOR_VLE, alpha=0.55, label="VLE")
        ax.set_xlim(low, high)
        if score.higher == "tie":
            panel_title = "%s\nAUC %.3f  no class difference" % (score.name, score.oriented)
        else:
            panel_title = "%s\nAUC %.3f  %s higher" % (score.name, score.oriented, score.higher)
        ax.set_title(panel_title, fontsize=7)
        ax.tick_params(labelsize=7)
        if panel % ncols == 0:
            ax.set_ylabel("Density", fontsize=8)
    for panel in range(n_features, flat.size):
        flat[panel].set_visible(False)
    fig.legend(
        handles=[
            Patch(facecolor=_COLOR_VLE, alpha=0.55, label="VLE"),
            Patch(facecolor=_COLOR_DECOY, alpha=0.55, label="DECOY"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.975),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Feature distributions by class, ordered by oriented AUC", fontsize=13, y=0.995)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    _save_fig(fig, path)


def _balance_weights(families: np.ndarray, is_vle: np.ndarray) -> np.ndarray:
    """Same family totals as construct_VLE_discriminator_score._row_weights."""
    py_dir = str(Path(__file__).resolve().parent)
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)
    from construct_VLE_discriminator_score import _row_weights

    return _row_weights(families, is_vle)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values, kind="mergesort")
    cumulative = np.cumsum(weights[order])
    total = float(cumulative[-1])
    if total <= 0.0:
        return float(values[order[0]])
    index = int(np.searchsorted(cumulative, quantile * total, side="left"))
    index = min(index, values.size - 1)
    return float(values[order[index]])


def _weighted_limits(column: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    low = _weighted_quantile(column, weights, 0.01)
    high = _weighted_quantile(column, weights, 0.99)
    if high > low:
        return low, high
    pad = 0.5 if low == 0.0 else max(0.5, abs(low) * 0.05)
    return low - pad, high + pad


def _weighted_auc(
    values: np.ndarray, is_vle: np.ndarray, weights: np.ndarray
) -> tuple[float, str]:
    """Weighted Mann-Whitney AUC. Above 0.5, VLE tends to have the larger value."""
    vle_mass = float(weights[is_vle].sum())
    decoy_mass = float(weights[~is_vle].sum())
    if vle_mass <= 0.0 or decoy_mass <= 0.0 or np.all(values == values[0]):
        return 0.5, "tie"
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    sorted_vle = is_vle[order]
    sorted_weights = weights[order]
    wins = 0.0
    decoy_below = 0.0
    start = 0
    n = values.size
    while start < n:
        end = start + 1
        while end < n and sorted_values[end] == sorted_values[start]:
            end += 1
        group_weights = sorted_weights[start:end]
        group_vle = sorted_vle[start:end]
        vle_weight = float(group_weights[group_vle].sum())
        decoy_weight = float(group_weights[~group_vle].sum())
        wins += vle_weight * decoy_below + 0.5 * vle_weight * decoy_weight
        decoy_below += decoy_weight
        start = end
    auc = wins / (vle_mass * decoy_mass)
    if auc > 0.5:
        return auc, "VLE"
    if auc < 0.5:
        return auc, "DECOY"
    return auc, "tie"


def _balanced_rows(
    feature_names: list[str],
    values: np.ndarray,
    is_vle: np.ndarray,
    weights: np.ndarray,
) -> list[tuple[str, float, float, str]]:
    rows = []
    for index, name in enumerate(feature_names):
        auc, higher = _weighted_auc(values[:, index], is_vle, weights)
        rows.append((name, auc, max(auc, 1.0 - auc), higher))
    rows.sort(key=lambda item: (-item[2], item[0]))
    return rows


def _balanced_auc_text(rows: list[tuple[str, float, float, str]]) -> str:
    lines = [
        "Family-balanced AUC of each circle feature.",
        "Each VLE family with at least 20 circles has the same total weight.",
        "VLE families below 20 keep weight 1 per row.",
        "Each DECOY family has the same total weight.",
        "Those totals sum to the VLE row count and the DECOY row count.",
        "AUC is P(weighted VLE > weighted DECOY) + 0.5 P(tie).",
        "Oriented AUC is max(AUC, 1-AUC). higher is the class with the larger values.",
        "",
        "feature\tauc\toriented_auc\thigher",
    ]
    for name, auc, oriented, higher in rows:
        lines.append("%s\t%.6g\t%.6g\t%s" % (name, auc, oriented, higher))
    lines.append("")
    return "\n".join(lines)


def _write_balanced_distributions(
    path: Path,
    rows: list[tuple[str, float, float, str]],
    feature_names: list[str],
    values: np.ndarray,
    is_vle: np.ndarray,
    weights: np.ndarray,
) -> None:
    n_features = len(rows)
    ncols = 6 if n_features > 6 else max(1, n_features)
    nrows = math.ceil(n_features / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.15 * ncols, 2.45 * nrows),
        squeeze=False,
    )
    flat = axes.ravel()
    name_to_index = {name: index for index, name in enumerate(feature_names)}
    for panel, (name, _auc, oriented, higher) in enumerate(rows):
        ax = flat[panel]
        column = values[:, name_to_index[name]]
        low, high = _weighted_limits(column, weights)
        bins = np.linspace(low, high, 41)
        keep = (column >= low) & (column <= high)
        decoy = keep & ~is_vle
        vle = keep & is_vle
        if np.any(decoy) and float(weights[decoy].sum()) > 0.0:
            ax.hist(
                column[decoy],
                bins=bins,
                weights=weights[decoy],
                density=True,
                color=_COLOR_DECOY,
                alpha=0.55,
                label="DECOY",
            )
        if np.any(vle) and float(weights[vle].sum()) > 0.0:
            ax.hist(
                column[vle],
                bins=bins,
                weights=weights[vle],
                density=True,
                color=_COLOR_VLE,
                alpha=0.55,
                label="VLE",
            )
        ax.set_xlim(low, high)
        if higher == "tie":
            panel_title = "%s\nAUC %.3f  no class difference" % (name, oriented)
        else:
            panel_title = "%s\nAUC %.3f  %s higher" % (name, oriented, higher)
        ax.set_title(panel_title, fontsize=7)
        ax.tick_params(labelsize=7)
        if panel % ncols == 0:
            ax.set_ylabel("Density", fontsize=8)
    for panel in range(n_features, flat.size):
        flat[panel].set_visible(False)
    fig.legend(
        handles=[
            Patch(facecolor=_COLOR_VLE, alpha=0.55, label="VLE"),
            Patch(facecolor=_COLOR_DECOY, alpha=0.55, label="DECOY"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.975),
        ncol=2,
        frameon=False,
    )
    fig.suptitle(
        "Feature distributions by class, equal family weight, ordered by oriented AUC",
        fontsize=13,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    _save_fig(fig, path)


def _spearman(values: np.ndarray) -> np.ndarray:
    n_features = values.shape[1]
    centered = np.empty_like(values)
    spread = np.empty(n_features, dtype=np.float64)
    for feature_index in range(n_features):
        ranks = _average_ranks(values[:, feature_index])
        centered[:, feature_index] = ranks - ranks.mean()
        spread[feature_index] = math.sqrt(float(np.dot(centered[:, feature_index], centered[:, feature_index])))
    corr = np.full((n_features, n_features), np.nan, dtype=np.float64)
    for i in range(n_features):
        if spread[i] == 0.0:
            continue
        for j in range(i, n_features):
            if spread[j] == 0.0:
                continue
            rho = float(np.dot(centered[:, i], centered[:, j]) / (spread[i] * spread[j]))
            corr[i, j] = rho
            corr[j, i] = rho
    np.fill_diagonal(corr, 1.0)
    return corr


def _write_correlation(
    path: Path,
    scores: list[_Score],
    feature_names: list[str],
    values: np.ndarray,
) -> None:
    name_to_index = {name: index for index, name in enumerate(feature_names)}
    order = [name_to_index[score.name] for score in scores]
    ordered = values[:, order]
    names = [score.name for score in scores]
    corr = _spearman(ordered)
    side = max(9.5, 0.36 * len(names) + 2.0)
    fig, ax = plt.subplots(figsize=(side + 2.0, side))
    cmap = plt.cm.coolwarm.copy()
    cmap.set_bad("0.85")
    image = ax.imshow(
        np.ma.masked_invalid(corr),
        vmin=-1.0,
        vmax=1.0,
        cmap=cmap,
        interpolation="nearest",
        aspect="equal",
    )
    positions = np.arange(len(names))
    ax.set_xticks(positions)
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_title("Spearman correlation, features ordered by oriented AUC")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Spearman rho")
    fig.tight_layout()
    _save_fig(fig, path)


def main() -> None:
    if len(sys.argv) > 3:
        _die("ERROR: usage: analyze_circle_features.py [input.tsv] [output_dir]")
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else _DEFAULT_INPUT
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else _DEFAULT_OUTPUT
    feature_names, values, is_vle, families = _load(input_path)
    scores = _score_features(feature_names, values, is_vle)
    adjusted = [score for score in scores if score.oriented_adj is not None]
    adjusted.sort(key=lambda score: (-(score.oriented_adj or 0.0), -abs(score.cohens_d), score.name))
    near_constant = _near_constant(feature_names, values)
    output_dir.mkdir(parents=True, exist_ok=True)
    ranking_path = output_dir / "feature_ranking.tsv"
    summary_path = output_dir / "summary.txt"
    _write_ranking(ranking_path, scores)
    text = _summary_text(int(is_vle.sum()), int((~is_vle).sum()), scores, adjusted, near_constant)
    _write_summary(summary_path, text)
    _write_ranking_bars(
        output_dir / "ranking_auc.png",
        scores,
        "Oriented AUC of each circle feature",
        adjusted=False,
    )
    _write_ranking_bars(
        output_dir / "ranking_auc_length_adjusted.png",
        scores,
        "Oriented AUC after residualizing on log(%s)" % _LENGTH_FEATURE,
        adjusted=True,
    )
    _write_distributions(
        output_dir / "distributions.png",
        scores,
        feature_names,
        values,
        is_vle,
    )
    _write_correlation(
        output_dir / "spearman_correlation.png",
        scores,
        feature_names,
        values,
    )
    balance = _balance_weights(families, is_vle)
    balanced_rows = _balanced_rows(feature_names, values, is_vle, balance)
    balanced_auc_path = output_dir / "feature_auc_balanced.txt"
    balanced_auc_path.write_text(_balanced_auc_text(balanced_rows), encoding="utf-8", newline="\n")
    balanced_dist_path = output_dir / "distributions_balanced.png"
    _write_balanced_distributions(
        balanced_dist_path,
        balanced_rows,
        feature_names,
        values,
        is_vle,
        balance,
    )
    print(text, end="")
    print("wrote %s" % ranking_path)
    print("wrote %s" % summary_path)
    print("wrote %s" % (output_dir / "ranking_auc.png"))
    print("wrote %s" % (output_dir / "ranking_auc_length_adjusted.png"))
    print("wrote %s" % (output_dir / "distributions.png"))
    print("wrote %s" % (output_dir / "spearman_correlation.png"))
    print("wrote %s" % balanced_auc_path)
    print("wrote %s" % balanced_dist_path)


if __name__ == "__main__":
    main()
