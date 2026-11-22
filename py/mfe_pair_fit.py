#!/usr/bin/env python3
"""Predict MFE mean and sd from a pairing-product fit of the RNAfold grid.

cmean and cstddev are fit by weighted least squares on every grid row.
Features are base frequencies (fT dropped because frequencies sum to 1),
the pairing products fC*fG, fA*fT, and fG*fT, and 1/L. Predicted Z uses
the same scaling as mfe_grid_interp:

    pred_mean = cmean_hat * L / 100
    pred_sd = cstddev_hat * L / 100
    Z = (MFE - pred_mean) / pred_sd

Existing measured scores are reused. If the median of predicted Z minus
measured Z stays beyond 0.5, the squares of the three pairing products are
added and that global fit is repeated once.

One polynomial still sits above the MFE surface between composition-grid
points. When the squared-product fit is also beyond 0.5, the same features
are fit locally for each query. The weight of a grid row is

    n * exp(-0.5 * (||f - f_query||^2 / sigma_f^2 + log(L/L_query)^2 / sigma_L^2))

sigma_f is half the 0.125 frequency step. sigma_L is one 1.2-fold length
step. Those bandwidths come from the grid spacing, not from the test scores.
"""

from __future__ import annotations

import math
import statistics
import sys
from pathlib import Path

import numpy as np

import correlate_measured_predicted_Z as corr
import permute_fold

_ROOT = Path(__file__).resolve().parent.parent
_BASE_FEATURES = (
    "intercept",
    "fA",
    "fC",
    "fG",
    "fC_fG",
    "fA_fT",
    "fG_fT",
    "inv_L",
)
_SQUARE_FEATURES = ("fC_fG_sq", "fA_fT_sq", "fG_fT_sq")
# Half of the tabulated frequency step, and one step of the 1.2-fold length grid.
_SIGMA_COMPOSITION = 0.0625
_SIGMA_LOG_L = math.log(1.2)
_LOCAL_WEIGHT_FLOOR = 1e-6
_SCORE_COLUMNS = (
    "set",
    "label",
    "L",
    "fA",
    "fC",
    "fG",
    "fT",
    "MFE",
    "measured_mean",
    "measured_sd",
    "measured_Z",
    "pred_mean",
    "pred_sd",
    "predicted_Z",
    "status",
)


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _features(length, fA, fC, fG, fT, squares: bool) -> np.ndarray:
    pairing_cg = fC * fG
    pairing_at = fA * fT
    pairing_gt = fG * fT
    columns = [
        np.ones_like(length, dtype=np.float64),
        fA,
        fC,
        fG,
        pairing_cg,
        pairing_at,
        pairing_gt,
        1.0 / length,
    ]
    if squares:
        columns.extend([pairing_cg ** 2, pairing_at ** 2, pairing_gt ** 2])
    return np.column_stack(columns)


def _fit_weighted(features: np.ndarray, target: np.ndarray, n_samples: np.ndarray) -> np.ndarray:
    weight = np.sqrt(n_samples)
    coef, _, _, _ = np.linalg.lstsq(features * weight[:, None], target * weight, rcond=None)
    return coef


def load_grid(path: Path) -> dict[str, np.ndarray]:
    try:
        data = np.loadtxt(path, comments="#", dtype=np.float64, ndmin=2)
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    except ValueError as exc:
        _die(f"ERROR: cannot parse '{path}': {exc}")
    if data.size == 0 or data.shape[1] != 10:
        _die(f"ERROR: expected 10 columns in {path}")
    return {
        "L": data[:, 0],
        "fA": data[:, 1],
        "fC": data[:, 2],
        "fG": data[:, 3],
        "fT": data[:, 4],
        "n": data[:, 5],
        "cmean": data[:, 8],
        "cstddev": data[:, 9],
    }


def fit_grid(grid: dict[str, np.ndarray], squares: bool) -> tuple[np.ndarray, np.ndarray]:
    features = _features(grid["L"], grid["fA"], grid["fC"], grid["fG"], grid["fT"], squares)
    cmean = _fit_weighted(features, grid["cmean"], grid["n"])
    cstddev = _fit_weighted(features, grid["cstddev"], grid["n"])
    return cmean, cstddev


def predict_row(
    length: float,
    fA: float,
    fC: float,
    fG: float,
    fT: float,
    mfe: float | None,
    cmean_coef: np.ndarray,
    csd_coef: np.ndarray,
    squares: bool,
) -> tuple[float, float, float | None]:
    point = _features(
        np.array([length]),
        np.array([fA]),
        np.array([fC]),
        np.array([fG]),
        np.array([fT]),
        squares,
    )
    cmean = (point @ cmean_coef).item()
    cstddev = (point @ csd_coef).item()
    pred_mean = cmean * length / 100.0
    pred_sd = cstddev * length / 100.0
    if mfe is None:
        return pred_mean, pred_sd, None
    return pred_mean, pred_sd, permute_fold.z_score(mfe, pred_mean, pred_sd)


def _parse_float(text: str) -> float | None:
    if text == "":
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def read_scores(path: Path, set_name: str) -> list[dict[str, object]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            index = {name: i for i, name in enumerate(header)}
            rows = []
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < len(header):
                    fields.extend([""] * (len(header) - len(fields)))
                status = fields[index["status"]]
                row: dict[str, object] = {
                    "set": set_name,
                    "label": fields[index["label"]],
                    "status": status,
                    "status_file": status,
                }
                row["L"] = _parse_float(fields[index["L"]])
                for name in ("fA", "fC", "fG", "fT", "MFE", "measured_mean", "measured_sd", "measured_Z", "pred_mean", "pred_sd", "predicted_Z"):
                    row[name] = _parse_float(fields[index[name]])
                row["old_predicted_Z"] = row["predicted_Z"]
                rows.append(row)
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    except KeyError as exc:
        _die(f"ERROR: {path} is missing column {exc}")
    return rows


def _clear_prediction(row: dict[str, object]) -> None:
    row["status"] = row["status_file"]
    row["pred_mean"] = None
    row["pred_sd"] = None
    row["predicted_Z"] = None


def _store_prediction(
    row: dict[str, object],
    pred_mean: float,
    pred_sd: float,
    predicted_z: float | None,
) -> None:
    row["status"] = row["status_file"]
    row["pred_mean"] = pred_mean
    row["pred_sd"] = pred_sd
    row["predicted_Z"] = predicted_z
    if isinstance(row["MFE"], float) and predicted_z is None:
        row["status"] = "zero_pred_sd"


def _query_frequencies(row: dict[str, object]) -> tuple[float, float, float, float] | None:
    length = row["L"]
    if not isinstance(length, float) or length <= 0.0:
        return None
    freqs = [row[name] for name in ("fA", "fC", "fG", "fT")]
    if any(not isinstance(value, float) for value in freqs):
        return None
    return float(freqs[0]), float(freqs[1]), float(freqs[2]), float(freqs[3])


def apply_fit(rows: list[dict[str, object]], cmean_coef: np.ndarray, csd_coef: np.ndarray, squares: bool) -> None:
    for row in rows:
        freqs = _query_frequencies(row)
        if freqs is None:
            _clear_prediction(row)
            continue
        mfe = row["MFE"] if isinstance(row["MFE"], float) else None
        pred_mean, pred_sd, predicted_z = predict_row(
            float(row["L"]),
            freqs[0],
            freqs[1],
            freqs[2],
            freqs[3],
            mfe,
            cmean_coef,
            csd_coef,
            squares,
        )
        _store_prediction(row, pred_mean, pred_sd, predicted_z)


def apply_local_fit(rows: list[dict[str, object]], grid: dict[str, np.ndarray], squares: bool) -> None:
    """Fit the pairing-product features separately around each query."""
    features = _features(grid["L"], grid["fA"], grid["fC"], grid["fG"], grid["fT"], squares)
    length_grid = grid["L"]
    frequency_grid = (grid["fA"], grid["fC"], grid["fG"], grid["fT"])
    n_samples = grid["n"]
    for row in rows:
        freqs = _query_frequencies(row)
        if freqs is None:
            _clear_prediction(row)
            continue
        length = float(row["L"])
        distance2 = sum((column - freq) ** 2 for column, freq in zip(frequency_grid, freqs))
        log_length = np.log(length_grid / length)
        raw_weight = n_samples * np.exp(
            -0.5
            * (
                distance2 / _SIGMA_COMPOSITION**2
                + log_length**2 / _SIGMA_LOG_L**2
            )
        )
        keep = raw_weight > raw_weight.max() * _LOCAL_WEIGHT_FLOOR
        weight = np.sqrt(raw_weight[keep])
        design = features[keep] * weight[:, None]
        cmean_coef = np.linalg.lstsq(design, grid["cmean"][keep] * weight, rcond=None)[0]
        csd_coef = np.linalg.lstsq(design, grid["cstddev"][keep] * weight, rcond=None)[0]
        mfe = row["MFE"] if isinstance(row["MFE"], float) else None
        pred_mean, pred_sd, predicted_z = predict_row(
            length,
            freqs[0],
            freqs[1],
            freqs[2],
            freqs[3],
            mfe,
            cmean_coef,
            csd_coef,
            squares,
        )
        _store_prediction(row, pred_mean, pred_sd, predicted_z)


def _finite_pair_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        row
        for row in rows
        if isinstance(row["measured_Z"], float)
        and isinstance(row["old_predicted_Z"], float)
        and isinstance(row["predicted_Z"], float)
    ]


def _bias(measured: list[float], predicted: list[float]) -> tuple[float, float]:
    if not measured:
        return float("nan"), float("nan")
    residual = [pred - meas for meas, pred in zip(measured, predicted)]
    right = sum(1 for meas, pred in zip(measured, predicted) if meas > pred)
    return statistics.median(residual), right / len(measured)


def _bias_of(rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], float, float, float, float]:
    shared = _finite_pair_rows(rows)
    measured = [float(row["measured_Z"]) for row in shared]
    old_median, old_right = _bias(measured, [float(row["old_predicted_Z"]) for row in shared])
    new_median, new_right = _bias(measured, [float(row["predicted_Z"]) for row in shared])
    return shared, old_median, old_right, new_median, new_right


def _fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.6g}"
    return str(value)


def _write_scores(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(_SCORE_COLUMNS) + "\n")
        for row in rows:
            length = row["L"]
            if isinstance(length, float) and length == round(length):
                row_text = dict(row)
                row_text["L"] = str(int(length))
            else:
                row_text = row
            handle.write("\t".join(_fmt(row_text[name]) for name in _SCORE_COLUMNS) + "\n")


def _write_coefficients(path: Path, names: tuple[str, ...], cmean: np.ndarray, cstddev: np.ndarray) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("feature\tcmean\tcstddev\n")
        for name, mean_coef, sd_coef in zip(names, cmean, cstddev):
            handle.write("%s\t%.6g\t%.6g\n" % (name, float(mean_coef), float(sd_coef)))


def _write_local_specification(path: Path, names: tuple[str, ...]) -> None:
    """Record bandwidths. Local coefficients depend on the query, so they are not stored."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("feature\tcmean\tcstddev\n")
        handle.write("sigma_composition\t%.6g\t%.6g\n" % (_SIGMA_COMPOSITION, _SIGMA_COMPOSITION))
        handle.write("sigma_log_L\t%.6g\t%.6g\n" % (_SIGMA_LOG_L, _SIGMA_LOG_L))
        for name in names:
            handle.write("%s\t\t\n" % name)


def _plot_rows(rows: list[dict[str, object]], zmin: float | None) -> tuple[tuple[list[float], list[float], list[str]], tuple[list[float], list[float], list[str]]]:
    controls: tuple[list[float], list[float], list[str]] = ([], [], [])
    rods: tuple[list[float], list[float], list[str]] = ([], [], [])
    for row in rows:
        measured_z = row["measured_Z"]
        predicted_z = row["predicted_Z"]
        if not isinstance(measured_z, float) or not isinstance(predicted_z, float):
            continue
        if zmin is not None and (measured_z < zmin or predicted_z < zmin):
            continue
        target = controls if row["set"] == "controls" else rods
        target[0].append(measured_z)
        target[1].append(predicted_z)
        target[2].append(str(row["label"]))
    return controls, rods


def _vle_axis_limits(controls, rods) -> tuple[float, float]:
    """Lowest VLE Z on either axis, through the highest Z in the plotted data."""
    measured, predicted, labels = controls
    vle_values = [
        value
        for x, y, label in zip(measured, predicted, labels)
        if str(label).startswith("VLE")
        for value in (x, y)
    ]
    if not vle_values:
        _die("ERROR: no VLE scores to set the scatter axis")
    all_values = list(measured) + list(predicted) + list(rods[0]) + list(rods[1])
    return min(vle_values), max(all_values)


def _save_plot(path: Path, controls, rods, limits: tuple[float, float] | None = None) -> tuple[float, float, float]:
    temporary = path.with_suffix(".tmp.png")
    stats = corr._write_combined(temporary, controls, rods, limits)
    temporary.replace(path)
    return stats


def main() -> None:
    grid_path = _ROOT / "output" / "rnafold_grid_mean_var.tsv"
    outdir = _ROOT / "correlate_measured_predicted_Z"
    controls = read_scores(outdir / "controls_z_scores.tsv", "controls")
    rods = read_scores(outdir / "ROD_z_scores.tsv", "RODs")
    rows = controls + rods
    grid = load_grid(grid_path)

    squares = False
    local = False
    feature_names = _BASE_FEATURES
    cmean_coef, csd_coef = fit_grid(grid, squares)
    apply_fit(rows, cmean_coef, csd_coef, squares)
    shared, old_median, old_right, new_median, new_right = _bias_of(rows)
    print(
        "global base median predicted_Z - measured_Z = %.6g  right = %.3f (n = %d)"
        % (new_median, new_right, len(shared))
    )
    if abs(new_median) > 0.5:
        squares = True
        feature_names = _BASE_FEATURES + _SQUARE_FEATURES
        cmean_coef, csd_coef = fit_grid(grid, squares)
        apply_fit(rows, cmean_coef, csd_coef, squares)
        shared, old_median, old_right, new_median, new_right = _bias_of(rows)
        print(
            "global squares median predicted_Z - measured_Z = %.6g  right = %.3f (n = %d)"
            % (new_median, new_right, len(shared))
        )
    if abs(new_median) > 0.5:
        local = True
        apply_local_fit(rows, grid, squares)
        shared, old_median, old_right, new_median, new_right = _bias_of(rows)
        print(
            "local fit sigma_composition = %.6g  sigma_log_L = %.6g"
            % (_SIGMA_COMPOSITION, _SIGMA_LOG_L)
        )

    if local:
        _write_local_specification(outdir / "pairfit_coefficients.tsv", feature_names)
    else:
        _write_coefficients(outdir / "pairfit_coefficients.tsv", feature_names, cmean_coef, csd_coef)
    _write_scores(outdir / "pairfit_z_scores.tsv", rows)
    full_controls, full_rods = _plot_rows(rows, None)
    clipped_controls, clipped_rods = _plot_rows(rows, -10.0)
    axis_limits = _vle_axis_limits(full_controls, full_rods)
    full_r = _save_plot(outdir / "scatter_pairfit.png", full_controls, full_rods, axis_limits)
    clipped_r = _save_plot(outdir / "scatter_pairfit_Z_-10_to_max.png", clipped_controls, clipped_rods)
    print("scatter axis %.6g to %.6g" % axis_limits)

    print("accepted %s features %s" % ("local" if local else "global", " ".join(feature_names)))
    print("interpolation median predicted_Z - measured_Z = %.6g  right = %.3f (n = %d)" % (old_median, old_right, len(shared)))
    print("accepted median predicted_Z - measured_Z = %.6g  right = %.3f (n = %d)" % (new_median, new_right, len(shared)))
    print(
        "all pairs Pearson r RODs = %.6g (n = %d) controls = %.6g (n = %d) both = %.6g (n = %d)"
        % (full_r[0], len(full_rods[0]), full_r[1], len(full_controls[0]), full_r[2], len(full_rods[0]) + len(full_controls[0]))
    )
    print(
        "Z>=-10 Pearson r RODs = %.6g (n = %d) controls = %.6g (n = %d) both = %.6g (n = %d)"
        % (
            clipped_r[0],
            len(clipped_rods[0]),
            clipped_r[1],
            len(clipped_controls[0]),
            clipped_r[2],
            len(clipped_rods[0]) + len(clipped_controls[0]),
        )
    )


if __name__ == "__main__":
    main()
