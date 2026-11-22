#!/usr/bin/env python3
"""Predict control Z-scores from length, composition, and RNAfold MFE.

The null mean and sd are the local pairing-product fit in mfe_pair_fit.
Z = (MFE - pred_mean) / pred_sd, with pred_mean = cmean_hat * L / 100.
"""

from __future__ import annotations

import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import mfe_pair_fit as fit

_ROOT = Path(__file__).resolve().parent.parent
_COLUMNS = (
    "label",
    "group",
    "class",
    "L",
    "NA",
    "NC",
    "NG",
    "NT",
    "MFE",
    "pred_mean",
    "pred_sd",
    "predicted_Z",
)
_N_BINS = 32
_CHUNK = 2500
_WORKERS = min(32, os.cpu_count() or 1)

_GRID: dict[str, np.ndarray] = {}


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _load_mfe(path: Path) -> dict[str, str]:
    scores: dict[str, str] = {}
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for line in handle:
                text = line.rstrip("\r\n")
                if not text:
                    continue
                label, mfe = text.split("\t")
                scores[label] = mfe
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    except ValueError:
        _die(f"ERROR: expected label and MFE in {path}")
    return scores


def _group_class(label: str) -> tuple[str, str]:
    group, _, rest = label.partition("_")
    class_name, _, _rest = rest.partition("_")
    return group, class_name


def _read_controls(path: Path, mfe_text: dict[str, str]) -> tuple[list[dict[str, object]], np.ndarray]:
    """Return row records and an (n, 5) array of L, fA, fC, fG, fT."""
    rows: list[dict[str, object]] = []
    queries: list[list[float]] = []
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for line_number, line in enumerate(handle, start=1):
                text = line.rstrip("\r\n")
                if not text:
                    continue
                fields = text.split("\t")
                if len(fields) != 6:
                    _die(f"ERROR: {path}:{line_number} expected 6 columns")
                label, length_text, na_text, nc_text, ng_text, nt_text = fields
                group, class_name = _group_class(label)
                if group not in ("VLE", "DECOY") or not class_name:
                    _die(f"ERROR: {path}:{line_number} label is not VLE or DECOY: {label}")
                try:
                    length = int(length_text)
                    counts = (int(na_text), int(nc_text), int(ng_text), int(nt_text))
                except ValueError:
                    _die(f"ERROR: {path}:{line_number} has a non-integer length or count")
                record: dict[str, object] = {
                    "label": label,
                    "group": group,
                    "class": class_name,
                    "L": length,
                    "NA": counts[0],
                    "NC": counts[1],
                    "NG": counts[2],
                    "NT": counts[3],
                    "MFE": mfe_text.get(label),
                }
                if length > 0:
                    scale = float(length)
                    queries.append(
                        [
                            scale,
                            counts[0] / scale,
                            counts[1] / scale,
                            counts[2] / scale,
                            counts[3] / scale,
                        ]
                    )
                    record["query"] = len(queries) - 1
                rows.append(record)
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    if not queries:
        _die(f"ERROR: no positive-length controls in {path}")
    return rows, np.asarray(queries, dtype=np.float64)


def _init_worker(grid_path: str) -> None:
    grid = fit.load_grid(Path(grid_path))
    _GRID["features"] = fit._features(grid["L"], grid["fA"], grid["fC"], grid["fG"], grid["fT"], True)
    _GRID["f"] = np.column_stack((grid["fA"], grid["fC"], grid["fG"], grid["fT"]))
    _GRID["log_L"] = np.log(grid["L"])
    _GRID["n"] = grid["n"]
    _GRID["y"] = np.column_stack((grid["cmean"], grid["cstddev"]))
    _GRID["sigma_f2"] = fit._SIGMA_COMPOSITION ** 2
    _GRID["sigma_l2"] = fit._SIGMA_LOG_L ** 2
    _GRID["floor"] = fit._LOCAL_WEIGHT_FLOOR


def _predict_chunk(queries: np.ndarray) -> np.ndarray:
    """Return (n, 2) predicted cmean and cstddev. NaN when the fit is empty."""
    features = _GRID["features"]
    frequency = _GRID["f"]
    log_length = _GRID["log_L"]
    n_samples = _GRID["n"]
    target = _GRID["y"]
    sigma_f2 = _GRID["sigma_f2"]
    sigma_l2 = _GRID["sigma_l2"]
    floor = _GRID["floor"]
    predicted = np.full((queries.shape[0], 2), np.nan, dtype=np.float64)
    for index in range(queries.shape[0]):
        length = queries[index, 0]
        freqs = queries[index, 1:]
        distance2 = np.sum((frequency - freqs) ** 2, axis=1)
        log_gap = log_length - math.log(length)
        raw_weight = n_samples * np.exp(-0.5 * (distance2 / sigma_f2 + (log_gap ** 2) / sigma_l2))
        limit = raw_weight.max() * floor
        keep = raw_weight > limit
        if not np.any(keep):
            continue
        weight = np.sqrt(raw_weight[keep])
        design = features[keep] * weight[:, None]
        rhs = target[keep] * weight[:, None]
        coef = np.linalg.lstsq(design, rhs, rcond=None)[0]
        fA, fC, fG, fT = freqs
        pairing_cg = fC * fG
        pairing_at = fA * fT
        pairing_gt = fG * fT
        point = np.array(
            [
                1.0,
                fA,
                fC,
                fG,
                pairing_cg,
                pairing_at,
                pairing_gt,
                1.0 / length,
                pairing_cg ** 2,
                pairing_at ** 2,
                pairing_gt ** 2,
            ]
        )
        predicted[index] = point @ coef
    return predicted


def _predict_all(queries: np.ndarray, grid_path: Path) -> np.ndarray:
    chunks = [queries[start : start + _CHUNK] for start in range(0, queries.shape[0], _CHUNK)]
    predicted = np.empty((queries.shape[0], 2), dtype=np.float64)
    done = 0
    with ProcessPoolExecutor(
        max_workers=_WORKERS,
        initializer=_init_worker,
        initargs=(str(grid_path),),
    ) as pool:
        for chunk, result in zip(chunks, pool.map(_predict_chunk, chunks)):
            predicted[done : done + chunk.shape[0]] = result
            done += chunk.shape[0]
            print(f"predicted {done}/{queries.shape[0]}", file=sys.stderr, flush=True)
    return predicted


def _fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.6g}"
    return str(value)


def _write_table(path: Path, rows: list[dict[str, object]], scaled: np.ndarray) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(_COLUMNS) + "\n")
        for row in rows:
            query = row.get("query")
            pred_mean = None
            pred_sd = None
            predicted_z = None
            if isinstance(query, int):
                cmean = float(scaled[query, 0])
                cstddev = float(scaled[query, 1])
                length = int(row["L"])
                if math.isfinite(cmean):
                    pred_mean = cmean * length / 100.0
                if math.isfinite(cstddev):
                    pred_sd = cstddev * length / 100.0
                mfe_text = row["MFE"]
                if (
                    isinstance(mfe_text, str)
                    and isinstance(pred_mean, float)
                    and isinstance(pred_sd, float)
                    and pred_sd > 0.0
                ):
                    predicted_z = (float(mfe_text) - pred_mean) / pred_sd
            fields = [
                row["label"],
                row["group"],
                row["class"],
                row["L"],
                row["NA"],
                row["NC"],
                row["NG"],
                row["NT"],
                row["MFE"],
                pred_mean,
                pred_sd,
                predicted_z,
            ]
            handle.write("\t".join(_fmt(field) for field in fields) + "\n")


def _collect_z(rows: list[dict[str, object]], scaled: np.ndarray) -> tuple[dict[str, list[float]], dict[str, str]]:
    grouped: dict[str, list[float]] = {"VLE": [], "DECOY": []}
    by_class: dict[str, list[float]] = {}
    class_group: dict[str, str] = {}
    for row in rows:
        query = row.get("query")
        mfe_text = row["MFE"]
        if not isinstance(query, int) or not isinstance(mfe_text, str):
            continue
        cstddev = float(scaled[query, 1])
        length = int(row["L"])
        pred_sd = cstddev * length / 100.0
        if not (math.isfinite(pred_sd) and pred_sd > 0.0):
            continue
        cmean = float(scaled[query, 0])
        if not math.isfinite(cmean):
            continue
        predicted_z = (float(mfe_text) - cmean * length / 100.0) / pred_sd
        if not math.isfinite(predicted_z):
            continue
        group = str(row["group"])
        class_name = str(row["class"])
        grouped[group].append(predicted_z)
        by_class.setdefault(class_name, []).append(predicted_z)
        class_group[class_name] = group
    panels: dict[str, list[float]] = {"All VLE": grouped["VLE"], "All Decoy": grouped["DECOY"]}
    for class_name in sorted(by_class):
        panels[class_name] = by_class[class_name]
    return panels, class_group


def _panels_from_values(
    grouped_values: list[tuple[str, str, float]],
) -> tuple[dict[str, list[float]], dict[str, str]]:
    """Build All VLE, All Decoy, then alphabetical class panels.

    Each item is (group, class, value).
    """
    grouped: dict[str, list[float]] = {"VLE": [], "DECOY": []}
    by_class: dict[str, list[float]] = {}
    class_group: dict[str, str] = {}
    for group, class_name, value in grouped_values:
        grouped.setdefault(group, []).append(value)
        by_class.setdefault(class_name, []).append(value)
        class_group[class_name] = group
    panels: dict[str, list[float]] = {"All VLE": grouped.get("VLE", []), "All Decoy": grouped.get("DECOY", [])}
    for class_name in sorted(by_class):
        panels[class_name] = by_class[class_name]
    return panels, class_group


def _collect_mfe_per_l(rows: list[dict[str, object]]) -> tuple[dict[str, list[float]], dict[str, str]]:
    values: list[tuple[str, str, float]] = []
    for row in rows:
        mfe_text = row["MFE"]
        length = row["L"]
        if not isinstance(mfe_text, str) or not isinstance(length, int) or length <= 0:
            continue
        ratio = float(mfe_text) / float(length)
        if not math.isfinite(ratio):
            continue
        values.append((str(row["group"]), str(row["class"]), ratio))
    return _panels_from_values(values)


def _collect_mfe_per_l_table(path: Path) -> tuple[dict[str, list[float]], dict[str, str]]:
    values: list[tuple[str, str, float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        group_col = header.index("group")
        class_col = header.index("class")
        length_col = header.index("L")
        mfe_col = header.index("MFE")
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if fields[mfe_col] == "" or fields[length_col] in ("", "0"):
                continue
            length = float(fields[length_col])
            if length <= 0.0:
                continue
            ratio = float(fields[mfe_col]) / length
            if math.isfinite(ratio):
                values.append((fields[group_col], fields[class_col], ratio))
    return _panels_from_values(values)


def _axis_limits(values: np.ndarray, full_range: bool) -> tuple[float, float]:
    if full_range:
        low = float(values.min())
        high = float(values.max())
    else:
        # A few extrapolated Z scores reach about -1000 and +1000. Binning
        # from the absolute minimum to maximum puts the real distribution in one bin.
        low = float(np.quantile(values, 0.0005))
        high = float(np.quantile(values, 0.9995))
    if high <= low:
        high = low + 1.0
    return low, high


def _write_histograms(
    path: Path,
    panels: dict[str, list[float]],
    class_group: dict[str, str],
    xlabel: str,
    full_range: bool,
) -> None:
    finite = np.asarray([value for values in panels.values() for value in values], dtype=np.float64)
    if finite.size == 0:
        _die("ERROR: no finite values to plot")
    low, high = _axis_limits(finite, full_range)
    edges = np.linspace(low, high, _N_BINS + 1)
    titles = list(panels)
    columns = 4
    rows = math.ceil(len(titles) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(14, 3.1 * rows), sharex=True)
    flat = np.atleast_1d(axes).ravel()
    for ax, title in zip(flat, titles):
        if title == "All VLE":
            color = "magenta"
        elif title == "All Decoy":
            color = "orange"
        else:
            color = "magenta" if class_group.get(title) == "VLE" else "orange"
        ax.hist(panels[title], bins=edges, color=color)
        ax.set_xlim(low, high)
        ax.set_title("%s (n = %d)" % (title, len(panels[title])))
        ax.set_ylabel("Count")
    for ax in flat[len(titles) :]:
        ax.set_visible(False)
    for ax in flat[max(0, len(flat) - columns) :]:
        if ax.get_visible():
            ax.set_xlabel(xlabel)
    fig.suptitle("Shared bins %.4g to %.4g" % (low, high), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    temporary = path.with_suffix(".tmp.png")
    fig.savefig(temporary, dpi=150)
    plt.close(fig)
    temporary.replace(path)


def main() -> None:
    outdir = _ROOT / "big_output"
    grid_path = _ROOT / "output" / "rnafold_grid_mean_var.tsv"
    mfe_text = _load_mfe(outdir / "controls_label_mfe.tsv")
    rows, queries = _read_controls(outdir / "controls.compos.tsv", mfe_text)
    print(
        "controls %d  with MFE %d  queries %d  workers %d"
        % (len(rows), sum(1 for row in rows if row["MFE"] is not None), queries.shape[0], _WORKERS),
        file=sys.stderr,
        flush=True,
    )
    scaled = _predict_all(queries, grid_path)
    table_path = outdir / "controls_predicted_z.tsv"
    _write_table(table_path, rows, scaled)
    panels, class_group = _collect_z(rows, scaled)
    figure_path = outdir / "controls_predicted_z_hist.png"
    _write_histograms(figure_path, panels, class_group, "Predicted Z", False)
    mfe_panels, mfe_groups = _collect_mfe_per_l(rows)
    mfe_figure = outdir / "controls_mfe_per_L_hist.png"
    _write_histograms(mfe_figure, mfe_panels, mfe_groups, "MFE/L", True)
    print("wrote %s" % mfe_figure)
    print("wrote %s" % table_path)
    print("wrote %s" % figure_path)
    for title, values in panels.items():
        print("%s\tn=%d" % (title, len(values)))


if __name__ == "__main__":
    main()
