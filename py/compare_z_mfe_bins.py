#!/usr/bin/env python3
"""Compare predicted Z and MFE/L as one-dimensional class features.

Each class has equal total weight. Contiguous bins are chosen on a training
half to maximize macro-averaged recall, then scored on the held-out half.
Both features are compared at the same number of bins.
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
_N_ATOMS = 256
_K_MAX = 16
_N_HOLDOUTS = 10
_HOLDOUT_SEED = 1


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _load(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    predicted_z: list[float] = []
    mfe_per_l: list[float] = []
    class_names_row: list[str] = []
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            class_col = header.index("class")
            length_col = header.index("L")
            mfe_col = header.index("MFE")
            z_col = header.index("predicted_Z")
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if fields[z_col] == "" or fields[mfe_col] == "" or fields[length_col] in ("", "0"):
                    continue
                length = float(fields[length_col])
                if not math.isfinite(length) or length <= 0.0:
                    continue
                z_value = float(fields[z_col])
                ratio = float(fields[mfe_col]) / length
                if not (math.isfinite(z_value) and math.isfinite(ratio)):
                    continue
                predicted_z.append(z_value)
                mfe_per_l.append(ratio)
                class_names_row.append(fields[class_col])
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    except ValueError as exc:
        _die(f"ERROR: {path} is missing a column or has a bad value: {exc}")
    if not predicted_z:
        _die(f"ERROR: no comparable rows in {path}")
    class_names = sorted(set(class_names_row))
    class_index = np.fromiter((class_names.index(name) for name in class_names_row), dtype=np.int32, count=len(class_names_row))
    return np.asarray(predicted_z, dtype=np.float64), np.asarray(mfe_per_l, dtype=np.float64), class_index, class_names


def _stratified_split(class_index: np.ndarray, n_classes: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    train_parts: list[np.ndarray] = []
    test_parts: list[np.ndarray] = []
    for class_id in range(n_classes):
        members = np.flatnonzero(class_index == class_id)
        if members.size < 2:
            _die(f"ERROR: class {class_id} has {members.size} rows, need at least 2")
        order = rng.permutation(members.size)
        members = members[order]
        n_test = members.size // 2
        test_parts.append(members[:n_test])
        train_parts.append(members[n_test:])
    return np.concatenate(train_parts), np.concatenate(test_parts)


def _atom_index(n_rows: int, n_atoms: int) -> np.ndarray:
    """Equal-count atom id for each rank from 0 to n_rows - 1."""
    n_atoms = min(n_atoms, n_rows)
    base, extra = divmod(n_rows, n_atoms)
    sizes = np.full(n_atoms, base, dtype=np.int32)
    sizes[:extra] += 1
    return np.repeat(np.arange(n_atoms, dtype=np.int32), sizes)


def _segment_reward(prefix: np.ndarray, n_class: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reward and winning class for every atom segment [start, end)."""
    right = prefix[:, None, :]
    left = prefix[:, :, None]
    weight = (right - left) / n_class[:, None, None]
    return np.max(weight, axis=0), np.argmax(weight, axis=0)


def _best_cuts(reward: np.ndarray, k_max: int) -> tuple[np.ndarray, np.ndarray]:
    """Return total reward and predecessor start for k bins covering atoms [0, end)."""
    n_bounds = reward.shape[0]
    n_atoms = n_bounds - 1
    k_max = min(k_max, n_atoms)
    total = np.full((k_max + 1, n_bounds), -np.inf)
    start = np.full((k_max + 1, n_bounds), -1, dtype=np.int32)
    total[1, 1:] = reward[0, 1:]
    start[1, 1:] = 0
    for k in range(2, k_max + 1):
        for end in range(k, n_atoms + 1):
            candidates = total[k - 1, k - 1 : end] + reward[k - 1 : end, end]
            choice = int(np.argmax(candidates))
            start[k, end] = k - 1 + choice
            total[k, end] = candidates[choice]
    return total, start


def _bin_of_atom(start: np.ndarray, k: int, n_atoms: int) -> np.ndarray:
    bounds = [n_atoms]
    end = n_atoms
    bins_left = k
    while bins_left > 0:
        begin = int(start[bins_left, end])
        bounds.append(begin)
        end = begin
        bins_left -= 1
    bounds.reverse()
    bin_id = np.empty(n_atoms, dtype=np.int32)
    for bin_index, (begin, end) in enumerate(zip(bounds[:-1], bounds[1:])):
        bin_id[begin:end] = bin_index
    return bin_id


def _assign_bins(train_values: np.ndarray, train_order: np.ndarray, atom_of_rank: np.ndarray, bin_of_atom: np.ndarray, values: np.ndarray) -> np.ndarray:
    sorted_values = train_values[train_order]
    rank = np.searchsorted(sorted_values, values, side="left")
    np.clip(rank, 0, train_values.size - 1, out=rank)
    return bin_of_atom[atom_of_rank[rank]]


def _macro_recall(assigned: np.ndarray, class_index: np.ndarray, winners: np.ndarray, n_classes: int) -> float:
    recalls = np.empty(n_classes, dtype=np.float64)
    for class_id in range(n_classes):
        members = assigned[class_index == class_id]
        recalls[class_id] = np.mean(winners[members] == class_id)
    return float(np.mean(recalls))


def _winners_from_bins(prefix: np.ndarray, n_class: np.ndarray, bin_of_atom: np.ndarray) -> np.ndarray:
    n_bins = int(bin_of_atom.max()) + 1
    atom_counts = prefix[:, 1:] - prefix[:, :-1]
    winners = np.empty(n_bins, dtype=np.int32)
    for bin_id in range(n_bins):
        counts = atom_counts[:, bin_of_atom == bin_id].sum(axis=1)
        winners[bin_id] = int(np.argmax(counts / n_class))
    return winners


def _fit_feature(values: np.ndarray, class_index: np.ndarray, n_classes: int, k_max: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit partitions for K = 1..k_max. Return reward totals, starts, order, atom ids, class counts."""
    n_class = np.bincount(class_index, minlength=n_classes).astype(np.float64)
    order = np.argsort(values, kind="mergesort")
    atom_of_rank = _atom_index(values.size, _N_ATOMS)
    n_atoms = int(atom_of_rank[-1]) + 1
    counts = np.zeros((n_classes, n_atoms), dtype=np.float64)
    np.add.at(counts, (class_index[order], atom_of_rank), 1.0)
    prefix = np.zeros((n_classes, n_atoms + 1), dtype=np.float64)
    np.cumsum(counts, axis=1, out=prefix[:, 1:])
    reward, _winner = _segment_reward(prefix, n_class)
    total, start = _best_cuts(reward, min(k_max, n_atoms))
    return total, start, order, atom_of_rank, n_class


def _score_holdout(
    values: np.ndarray,
    class_index: np.ndarray,
    train_index: np.ndarray,
    test_index: np.ndarray,
    n_classes: int,
) -> np.ndarray:
    train_values = values[train_index]
    train_class = class_index[train_index]
    total, start, order, atom_of_rank, n_class = _fit_feature(train_values, train_class, n_classes, _K_MAX)
    n_atoms = int(atom_of_rank[-1]) + 1
    k_limit = min(_K_MAX, n_atoms)
    scores = np.full(_K_MAX, np.nan)
    prefix = np.zeros((n_classes, n_atoms + 1), dtype=np.float64)
    counts = np.zeros((n_classes, n_atoms), dtype=np.float64)
    np.add.at(counts, (train_class[order], atom_of_rank), 1.0)
    np.cumsum(counts, axis=1, out=prefix[:, 1:])
    for k in range(1, k_limit + 1):
        bin_of_atom = _bin_of_atom(start, k, n_atoms)
        winners = _winners_from_bins(prefix, n_class, bin_of_atom)
        assigned = _assign_bins(train_values, order, atom_of_rank, bin_of_atom, values[test_index])
        scores[k - 1] = _macro_recall(assigned, class_index[test_index], winners, n_classes)
    return scores


def _holdout_scores(values: np.ndarray, class_index: np.ndarray, n_classes: int) -> np.ndarray:
    """Return (n_holdouts, K_MAX) test scores."""
    scores = np.zeros((_N_HOLDOUTS, _K_MAX), dtype=np.float64)
    for fold in range(_N_HOLDOUTS):
        rng = np.random.default_rng(_HOLDOUT_SEED + fold)
        train_index, test_index = _stratified_split(class_index, n_classes, rng)
        scores[fold] = _score_holdout(values, class_index, train_index, test_index, n_classes)
        print(
            "fold %d  K=1 %.4f  best %.4f at K=%d" % (fold + 1, scores[fold, 0], np.nanmax(scores[fold]), int(np.nanargmax(scores[fold])) + 1),
            file=sys.stderr,
            flush=True,
        )
    return scores


def _best_k(mean_score: np.ndarray) -> int:
    finite = np.where(np.isfinite(mean_score))[0]
    best = finite[np.argmax(mean_score[finite])]
    return int(best) + 1


def _partition_all(values: np.ndarray, class_index: np.ndarray, n_classes: int, k: int) -> np.ndarray:
    total, start, order, atom_of_rank, _n_class = _fit_feature(values, class_index, n_classes, k)
    del total
    n_atoms = int(atom_of_rank[-1]) + 1
    bin_of_atom = _bin_of_atom(start, k, n_atoms)
    assigned_ranked = bin_of_atom[atom_of_rank]
    assigned = np.empty(values.size, dtype=np.int32)
    assigned[order] = assigned_ranked
    return assigned


def _write_scores(path: Path, feature_names: list[str], scores: list[np.ndarray]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("feature\tK\tmean_test_score\tstdev_test_score\n")
        for name, fold_scores in zip(feature_names, scores):
            for k in range(1, _K_MAX + 1):
                column = fold_scores[:, k - 1]
                handle.write(
                    "%s\t%d\t%.6g\t%.6g\n" % (name, k, float(np.mean(column)), float(np.std(column, ddof=1)))
                )


def _stack_fractions(assigned: np.ndarray, class_index: np.ndarray, n_class: np.ndarray) -> np.ndarray:
    n_bins = int(assigned.max()) + 1
    n_classes = n_class.size
    fractions = np.zeros((n_bins, n_classes), dtype=np.float64)
    for class_id in range(n_classes):
        counts = np.bincount(assigned[class_index == class_id], minlength=n_bins)
        fractions[:, class_id] = counts / n_class[class_id]
    return fractions


def _bin_edges(values: np.ndarray, assigned: np.ndarray) -> np.ndarray:
    n_bins = int(assigned.max()) + 1
    lows = np.empty(n_bins, dtype=np.float64)
    highs = np.empty(n_bins, dtype=np.float64)
    for bin_id in range(n_bins):
        members = values[assigned == bin_id]
        lows[bin_id] = members.min()
        highs[bin_id] = members.max()
    edges = np.empty(n_bins + 1, dtype=np.float64)
    edges[0] = lows[0]
    for bin_id in range(n_bins - 1):
        if highs[bin_id] < lows[bin_id + 1]:
            edges[bin_id + 1] = 0.5 * (highs[bin_id] + lows[bin_id + 1])
        else:
            edges[bin_id + 1] = highs[bin_id]
    edges[-1] = highs[-1]
    span = edges[-1] - edges[0]
    nudge = span * 0.002 if span > 0.0 else 1.0
    for bin_id in range(n_bins):
        if edges[bin_id + 1] <= edges[bin_id]:
            edges[bin_id + 1] = edges[bin_id] + nudge
    return edges


def _draw_partition(
    ax: plt.Axes,
    values: np.ndarray,
    assigned: np.ndarray,
    class_index: np.ndarray,
    n_class: np.ndarray,
    colors: list,
    title: str,
    xlabel: str,
) -> list:
    fractions = _stack_fractions(assigned, class_index, n_class)
    edges = _bin_edges(values, assigned)
    handles = []
    n_classes = n_class.size
    for bin_id in range(fractions.shape[0]):
        order = np.argsort(-fractions[bin_id])
        bottom = 0.0
        width = edges[bin_id + 1] - edges[bin_id]
        for class_id in order:
            height = fractions[bin_id, class_id]
            if height <= 0.0:
                continue
            bar = ax.bar(
                edges[bin_id],
                height,
                width=width,
                bottom=bottom,
                align="edge",
                color=colors[class_id],
                linewidth=0,
            )
            bottom += height
            if bin_id == 0:
                handles.append(bar)
            elif len(handles) <= class_id:
                handles.append(bar)
    # Rebuild legend handles in class order, one per class that appears.
    legend_handles = []
    for class_id in range(n_classes):
        patch = plt.Rectangle((0, 0), 1, 1, color=colors[class_id])
        legend_handles.append(patch)
    view = _view_limits(values)
    if view is None:
        ax.set_xlim(edges[0], edges[-1])
    else:
        ax.set_xlim(view[0], view[1])
    ax.set_ylabel("Fraction of class")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    return legend_handles


def _view_limits(values: np.ndarray) -> tuple[float, float] | None:
    """Clip an axis when a few extreme scores would hide the other bins."""
    low, high = np.quantile(values, [0.001, 0.999])
    interior = float(high - low)
    span = float(values.max() - values.min())
    if interior <= 0.0 or span / interior < 8.0:
        return None
    return float(low), float(high)


def _draw_curve(ax: plt.Axes, feature_names: list[str], scores: list[np.ndarray], n_classes: int) -> None:
    ks = np.arange(1, _K_MAX + 1)
    for name, fold_scores in zip(feature_names, scores):
        mean_score = fold_scores.mean(axis=0)
        spread = fold_scores.std(axis=0, ddof=1)
        ax.errorbar(ks, mean_score, yerr=spread, marker="o", linewidth=1.5, capsize=3, label=name)
    ax.axhline(1.0 / n_classes, color="0.4", linestyle="--", linewidth=1, label="one bin (1/17)")
    ax.set_xlabel("Number of bins K")
    ax.set_ylabel("Mean holdout purity")
    ax.set_xticks(ks)
    ax.legend(loc="upper left")


def _write_figure(
    path: Path,
    features: list[tuple[str, np.ndarray]],
    class_index: np.ndarray,
    class_names: list[str],
    scores: list[np.ndarray],
    chosen_k: list[int],
) -> None:
    n_classes = len(class_names)
    n_class = np.bincount(class_index, minlength=n_classes).astype(np.float64)
    colors = [plt.get_cmap("tab20")(index) for index in range(n_classes)]
    fig = plt.figure(figsize=(14, 10))
    grid = fig.add_gridspec(2, 2, height_ratios=[2.2, 1.0])
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    curve = fig.add_subplot(grid[1, :])
    legend_handles = None
    for ax, (name, values), fold_scores, k in zip(axes, features, scores, chosen_k):
        assigned = _partition_all(values, class_index, n_classes, k)
        mean_score = float(fold_scores[:, k - 1].mean())
        legend_handles = _draw_partition(
            ax,
            values,
            assigned,
            class_index,
            n_class,
            colors,
            "%s    K = %d    holdout purity = %.3f" % (name, k, mean_score),
            name,
        )
    _draw_curve(curve, [name for name, _values in features], scores, n_classes)
    fig.legend(
        legend_handles,
        class_names,
        loc="center left",
        bbox_to_anchor=(1.0, 0.62),
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 0.86, 1))
    temporary = path.with_suffix(".tmp.png")
    fig.savefig(temporary, dpi=150, bbox_inches="tight")
    plt.close(fig)
    temporary.replace(path)


def main() -> None:
    table_path = _ROOT / "big_output" / "controls_predicted_z.tsv"
    predicted_z, mfe_per_l, class_index, class_names = _load(table_path)
    n_classes = len(class_names)
    print(
        "rows %d  classes %d" % (predicted_z.size, n_classes),
        file=sys.stderr,
        flush=True,
    )
    features = [("predicted Z", predicted_z), ("MFE/L", mfe_per_l)]
    scores = []
    for name, values in features:
        print(name, file=sys.stderr, flush=True)
        scores.append(_holdout_scores(values, class_index, n_classes))
    chosen_k = [_best_k(fold_scores.mean(axis=0)) for fold_scores in scores]
    outdir = _ROOT / "big_output"
    score_path = outdir / "z_vs_mfe_purity.tsv"
    _write_scores(score_path, [name for name, _values in features], scores)
    figure_path = outdir / "z_vs_mfe_purity.png"
    _write_figure(figure_path, features, class_index, class_names, scores, chosen_k)
    print("wrote %s" % score_path)
    print("wrote %s" % figure_path)
    for (name, _values), fold_scores, k in zip(features, scores, chosen_k):
        print(
            "%s best K = %d  mean = %.4f  stdev = %.4f"
            % (name, k, float(fold_scores[:, k - 1].mean()), float(fold_scores[:, k - 1].std(ddof=1)))
        )


if __name__ == "__main__":
    main()
