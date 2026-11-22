#!/usr/bin/env python3
"""Predict Z from a headerless 7-column MFE table.

Columns are label, length, MFE, NA, NC, NG, NT. The null is the local
pairing-product fit in mfe_pair_fit, with the squared pairing products:

    pred_mean = cmean_hat * L / 100
    pred_sd = cstddev_hat * L / 100
    Z = (MFE - pred_mean) / pred_sd

Output is label, MFE, length, Z. Label, MFE, and length are copied from
the input. Z is blank when length is not positive, MFE is not finite, or
pred_sd is not positive.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import mfe_pair_fit as fit

_ROOT = Path(__file__).resolve().parent.parent
_N_FIELDS = 7


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _read_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for line_number, line in enumerate(handle, start=1):
                text = line.rstrip("\r\n")
                if not text:
                    continue
                fields = text.split("\t")
                if len(fields) != _N_FIELDS:
                    _die(f"ERROR: {path}:{line_number} expected {_N_FIELDS} columns")
                label, length_text, mfe_text, na_text, nc_text, ng_text, nt_text = fields
                try:
                    length = int(length_text)
                    counts = (int(na_text), int(nc_text), int(ng_text), int(nt_text))
                    mfe = float(mfe_text)
                except ValueError:
                    _die(f"ERROR: {path}:{line_number} has a non-numeric length, MFE, or count")
                row: dict[str, object] = {
                    "label_text": label,
                    "mfe_text": mfe_text,
                    "length_text": length_text,
                    "L": float(length),
                    "MFE": mfe if math.isfinite(mfe) else None,
                    "status_file": "",
                }
                if length > 0:
                    scale = float(length)
                    row["fA"] = counts[0] / scale
                    row["fC"] = counts[1] / scale
                    row["fG"] = counts[2] / scale
                    row["fT"] = counts[3] / scale
                rows.append(row)
    except OSError as exc:
        _die(f"ERROR: cannot read '{path}': {exc.strerror or exc}")
    if not rows:
        _die(f"ERROR: no rows in {path}")
    return rows


def _z_text(row: dict[str, object]) -> str:
    predicted = row.get("predicted_Z")
    if not isinstance(predicted, float) or not math.isfinite(predicted):
        return ""
    return f"{predicted:.6g}"


def _write_table(path: Path, rows: list[dict[str, object]]) -> int:
    scored = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("label\tMFE\tlength\tZ\n")
        for row in rows:
            z_text = _z_text(row)
            if z_text:
                scored += 1
            handle.write(
                f"{row['label_text']}\t{row['mfe_text']}\t{row['length_text']}\t{z_text}\n"
            )
    return scored


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Predict Z from a 7-column MFE table.")
    parser.add_argument("--input", required=True, help="Headerless TSV: label length MFE NA NC NG NT")
    parser.add_argument("--output", required=True, help="TSV: label MFE length Z")
    parser.add_argument("--grid", help="RNAfold mean/variance grid (default: output/rnafold_grid_mean_var.tsv)")
    args = parser.parse_args(argv)
    grid_path = Path(args.grid) if args.grid else _ROOT / "output" / "rnafold_grid_mean_var.tsv"
    rows = _read_rows(Path(args.input))
    grid = fit.load_grid(grid_path)
    print(f"rows {len(rows)}", file=sys.stderr, flush=True)
    fit.apply_local_fit(rows, grid, True)
    output_path = Path(args.output)
    scored = _write_table(output_path, rows)
    print(f"scored {scored}", file=sys.stderr, flush=True)
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
