#!/usr/bin/env python3
"""Assign sequences to MFE-binning leaves and write Z scores.

Z = (MFE - mean_train) / sqrt(var_train), using the leaf the sequence falls
into. mean_train and var_train are the training statistics stored in the tree.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import mfe_binning_features as feat

HEADER_TOKENS = {"label", "length", "L", "MFE", "NA", "NC", "NG", "NT"}
_OUTPUT_HEADER = "label\tMFE\tmean\tvariance\tZ"


def _first_fields(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        line = handle.readline()
    if line.startswith("\ufeff"):
        line = line[1:]
    return [field.strip() for field in line.rstrip("\r\n").split("\t")]


def _column_index(fields: list[str] | None) -> dict[str, int]:
    """Map the columns we need to indexes. No header means the fixed field order."""
    if fields is None:
        return {"label": 0, "L": 1, "MFE": 2, "NA": 3, "NC": 4, "NG": 5, "NT": 6}
    index = {name: i for i, name in enumerate(fields)}

    def need(*names: str) -> int:
        for name in names:
            if name in index:
                return index[name]
        raise ValueError(f"missing column {' or '.join(names)}")

    return {
        "label": need("label"),
        "L": need("length", "L"),
        "MFE": need("MFE"),
        "NA": need("NA"),
        "NC": need("NC"),
        "NG": need("NG"),
        "NT": need("NT"),
    }


def read_sequences(path: str):
    """Read label, length, MFE, NA, NC, NG, NT. A header is optional."""
    fields = _first_fields(path)
    has_header = any(field in HEADER_TOKENS for field in fields)
    columns = _column_index(fields if has_header else None)
    skiprows = 1 if has_header else 0
    try:
        labels = np.loadtxt(
            path,
            delimiter="\t",
            usecols=(columns["label"],),
            skiprows=skiprows,
            dtype=str,
            comments=None,
            ndmin=1,
            encoding="utf-8",
        )
        numeric = np.loadtxt(
            path,
            delimiter="\t",
            usecols=(
                columns["L"],
                columns["MFE"],
                columns["NA"],
                columns["NC"],
                columns["NG"],
                columns["NT"],
            ),
            skiprows=skiprows,
            dtype=np.float64,
            comments=None,
            ndmin=2,
            encoding="utf-8",
        )
    except ValueError as exc:
        raise ValueError(f"could not read {path}: {exc}") from exc
    labels = np.atleast_1d(labels)
    if numeric.ndim == 1:
        numeric = numeric.reshape(1, -1)
    if labels.shape[0] != numeric.shape[0]:
        raise ValueError(f"{path} label column and numeric columns have different lengths")
    return (
        labels,
        numeric[:, 0],
        numeric[:, 1],
        numeric[:, 2],
        numeric[:, 3],
        numeric[:, 4],
        numeric[:, 5],
    )


def leaf_mean_var_z(tree: feat.Tree, length, mfe, na, nc, ng, nt):
    """Return ``(mean, variance, Z, n_undefined)`` aligned with the input rows.

    ``n_undefined`` counts rows with ``length <= 0`` or a non-finite MFE.
    Those rows, and rows whose leaf variance is not positive, have Z ``nan``.
    """
    length = np.asarray(length, dtype=np.float64).reshape(-1)
    mfe = np.asarray(mfe, dtype=np.float64).reshape(-1)
    na = np.asarray(na, dtype=np.float64).reshape(-1)
    nc = np.asarray(nc, dtype=np.float64).reshape(-1)
    ng = np.asarray(ng, dtype=np.float64).reshape(-1)
    nt = np.asarray(nt, dtype=np.float64).reshape(-1)
    n = int(length.shape[0])
    mean = np.full(n, np.nan, dtype=np.float64)
    var = np.full(n, np.nan, dtype=np.float64)
    z = np.full(n, np.nan, dtype=np.float64)
    ok = (length > 0.0) & np.isfinite(mfe)
    n_undefined = int(n - np.count_nonzero(ok))
    if not np.any(ok):
        return mean, var, z, n_undefined
    leaf = feat.assign(tree, feat.features(length[ok], na[ok], nc[ok], ng[ok], nt[ok]))
    mean_ok = tree.mean_train[leaf]
    var_ok = tree.var_train[leaf]
    mean[ok] = mean_ok
    var[ok] = var_ok
    scored = np.isfinite(mean_ok) & np.isfinite(var_ok) & (var_ok > 0.0)
    z_ok = np.full(mean_ok.shape, np.nan, dtype=np.float64)
    z_ok[scored] = (mfe[ok][scored] - mean_ok[scored]) / np.sqrt(var_ok[scored])
    z[ok] = z_ok
    return mean, var, z, n_undefined


def _fmt(value: float) -> str:
    if not np.isfinite(value):
        return "nan"
    return "%.3g" % value


def write_z(path: str, labels, mfe, mean, var, z) -> None:
    """Write ``label MFE mean variance Z``. ``path`` of ``-`` writes to stdout."""
    labels = np.asarray(labels)
    mfe = np.asarray(mfe, dtype=np.float64).reshape(-1)
    mean = np.asarray(mean, dtype=np.float64).reshape(-1)
    var = np.asarray(var, dtype=np.float64).reshape(-1)
    z = np.asarray(z, dtype=np.float64).reshape(-1)
    n = int(labels.shape[0])
    own_handle = path != "-"
    handle = open(path, "w", encoding="utf-8", newline="\n") if own_handle else sys.stdout
    try:
        handle.write(_OUTPUT_HEADER + "\n")
        for i in range(n):
            handle.write(
                f"{labels[i]}\t{_fmt(mfe[i])}\t{_fmt(mean[i])}\t{_fmt(var[i])}\t{_fmt(z[i])}\n"
            )
    finally:
        if own_handle:
            handle.close()


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Write per-sequence MFE Z scores from a binning tree.")
    parser.add_argument("--tree", required=True, help="Tree TSV from mfe_binning_build_bins.py")
    parser.add_argument("--input", required=True, help="TSV: label length MFE NA NC NG NT, header optional")
    parser.add_argument("--output", default="-", help="Output TSV, or - for stdout")
    args = parser.parse_args(argv)
    try:
        tree = feat.load_tree(args.tree)
        labels, length, mfe, na, nc, ng, nt = read_sequences(args.input)
        mean, var, z, n_undefined = leaf_mean_var_z(tree, length, mfe, na, nc, ng, nt)
    except ValueError as exc:
        _die(str(exc))
    if n_undefined:
        print(
            f"undefined={n_undefined} (length<=0 or MFE not finite)",
            file=sys.stderr,
        )
    write_z(args.output, labels, mfe, mean, var, z)


if __name__ == "__main__":
    main()
