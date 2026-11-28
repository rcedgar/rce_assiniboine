#!/usr/bin/env python3
"""Circle-feature distributions for VLE_avsunviroid versus the selected decoys.

Same panels as analyze_circle_features.py distributions.png: density
histograms, features ordered by oriented Mann-Whitney AUC. The avsun class
is every row of big_avsum_thumbs/misc/avsun.tsv. The decoy class is every
row of the length-matched, per-family-capped decoy.tsv, including other
VLE families that were selected as decoys.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

from analyze_circle_features import (
    _COLOR_DECOY,
    _COLOR_VLE,
    _panel_limits,
    _score_features,
)
from dotbracket_circular_features import extract_circular_rna_features

_ROOT = Path(__file__).resolve().parent.parent
_AVSUN = _ROOT / "big_avsum_thumbs" / "misc" / "avsun.tsv"
_DECOY = _ROOT / "big_avsum_thumbs" / "misc" / "decoy.tsv"
_THUMBS = _ROOT / "big_avsum_thumbs" / "avsun_thumbnails"
_RANGE = _ROOT / "big_avsum_thumbs" / "misc" / "length_range.txt"
_OUT = _ROOT / "big_avsum_thumbs" / "distributions.png"
_COLOR_AVSUN = _COLOR_VLE
_COLOR_SELECTED = _COLOR_DECOY
_STEM_KEEP = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._")


def _die(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def _discard_avsun(feats: dict[str, float]) -> tuple[bool, bool, bool]:
    """True when pct_bulge > 10, pH < 17, or bulge_to_stem_ratio > 0.14."""
    bulge = feats["pct_bulge"] > 10
    ph = feats["pH"] < 17
    ratio = feats["bulge_to_stem_ratio"] > 0.14
    return bulge, ph, ratio


def _stem(label: str, line_no: int) -> str:
    stem = "".join(char if char in _STEM_KEEP else "_" for char in label)
    if len(stem) > 80:
        stem = stem[:80]
    folded = stem.upper()
    device = folded in ("CON", "PRN", "AUX", "NUL") or (
        len(folded) == 4 and folded[:3] in ("COM", "LPT") and folded[3] in "123456789"
    )
    if stem in ("", ".", "..") or device:
        return "record_%d" % line_no
    return stem


def _note_range(kept: int, discarded: int) -> None:
    if not _RANGE.is_file():
        return
    lines = _RANGE.read_text(encoding="utf-8").splitlines()
    out: list[str] = []
    saw_discarded = False
    for line in lines:
        if line.startswith("avsun_discarded "):
            out.append("avsun_discarded %d" % discarded)
            saw_discarded = True
        elif line.startswith("avsun "):
            out.append("avsun %d" % kept)
        else:
            out.append(line)
    if not saw_discarded:
        out.append("avsun_discarded %d" % discarded)
    _RANGE.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


def _filter_avsun() -> tuple[int, int]:
    """Rewrite avsun.tsv without discarded rows and delete their thumbnails."""
    kept_lines: list[str] = []
    kept_stems: set[str] = set()
    discarded_n = 0
    bulge_n = 0
    ph_n = 0
    ratio_n = 0
    with _AVSUN.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if line == "":
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                _die("ERROR: %s:%d expected label, sequence, and structure" % (_AVSUN, line_number))
            label, sequence, structure = fields[0], fields[1], fields[2]
            feats = extract_circular_rna_features(
                structure,
                max_unpaired_gap=5,
                train_median_bulge_len=1.0,
                sequence=sequence,
            )
            bulge, ph, ratio = _discard_avsun(feats)
            if bulge or ph or ratio:
                discarded_n += 1
                bulge_n += int(bulge)
                ph_n += int(ph)
                ratio_n += int(ratio)
                continue
            kept_lines.append(line)
            kept_stems.add(_stem(label, line_number))
            if line_number % 2000 == 0:
                print("filter avsun %d" % line_number, file=sys.stderr)
    if not kept_lines:
        _die("ERROR: every avsun row was discarded")
    with _AVSUN.open("w", encoding="utf-8", newline="\n") as handle:
        for line in kept_lines:
            handle.write(line + "\n")
    removed_png = 0
    if _THUMBS.is_dir():
        for png in _THUMBS.glob("*.png"):
            if png.name.startswith("."):
                continue
            if png.name not in {stem + ".png" for stem in kept_stems}:
                png.unlink()
                removed_png += 1
    # A row can fail more than one test, so the reason counts can overlap.
    print(
        "discarded %d avsun (pct_bulge %d, pH %d, bulge_to_stem_ratio %d), kept %d, removed %d png"
        % (discarded_n, bulge_n, ph_n, ratio_n, len(kept_lines), removed_png),
        file=sys.stderr,
    )
    _note_range(len(kept_lines), discarded_n)
    return len(kept_lines), discarded_n


def _load(path: Path, feature_names: list[str]) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if line == "":
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                _die("ERROR: %s:%d expected label, sequence, and structure" % (path, line_number))
            label, sequence, structure = fields[0], fields[1], fields[2]
            if len(sequence) != len(structure):
                _die("ERROR: %s:%d sequence and structure lengths differ: %s" % (path, line_number, label))
            feats = extract_circular_rna_features(
                structure,
                max_unpaired_gap=5,
                train_median_bulge_len=1.0,
                sequence=sequence,
            )
            rows.append([float(feats[name]) for name in feature_names])
            if line_number % 2000 == 0:
                print("features %s %d" % (path.name, line_number), file=sys.stderr)
    if not rows:
        _die("ERROR: no rows in %s" % path)
    print("features %s %d" % (path.name, len(rows)), file=sys.stderr)
    return np.asarray(rows, dtype=np.float64)


def _class_name(higher: str) -> str:
    if higher == "VLE":
        return "avsun"
    if higher == "DECOY":
        return "decoy"
    return higher


def _write_plot(path: Path, feature_names: list[str], values: np.ndarray, is_avsun: np.ndarray) -> None:
    scores = _score_features(feature_names, values, is_avsun)
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
        decoy = column[~is_avsun]
        avsun = column[is_avsun]
        decoy = decoy[(decoy >= low) & (decoy <= high)]
        avsun = avsun[(avsun >= low) & (avsun <= high)]
        if decoy.size:
            ax.hist(decoy, bins=bins, density=True, color=_COLOR_SELECTED, alpha=0.55, label="decoy")
        if avsun.size:
            ax.hist(avsun, bins=bins, density=True, color=_COLOR_AVSUN, alpha=0.55, label="avsun")
        ax.set_xlim(low, high)
        higher = _class_name(score.higher)
        if higher == "tie":
            panel_title = "%s\nAUC %.3f  no class difference" % (score.name, score.oriented)
        else:
            panel_title = "%s\nAUC %.3f  %s higher" % (score.name, score.oriented, higher)
        ax.set_title(panel_title, fontsize=7)
        ax.tick_params(labelsize=7)
        if panel % ncols == 0:
            ax.set_ylabel("Density", fontsize=8)
    for panel in range(n_features, flat.size):
        flat[panel].set_visible(False)
    fig.legend(
        handles=[
            Patch(facecolor=_COLOR_AVSUN, alpha=0.55, label="avsun"),
            Patch(facecolor=_COLOR_SELECTED, alpha=0.55, label="decoy"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.975),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Feature distributions by class, ordered by oriented AUC", fontsize=13, y=0.995)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    if not _AVSUN.is_file():
        _die("ERROR: missing %s" % _AVSUN)
    filter_only = len(sys.argv) == 2 and sys.argv[1] == "--filter-only"
    if len(sys.argv) > 2 or (len(sys.argv) == 2 and not filter_only):
        _die("ERROR: usage: avsun_feature_distributions.py [--filter-only]")
    kept, discarded = _filter_avsun()
    if filter_only:
        print("%d %d" % (kept, discarded))
        return
    if not _DECOY.is_file():
        _die("ERROR: missing %s" % _DECOY)
    feature_names = sorted(extract_circular_rna_features("").keys())
    if not feature_names:
        _die("ERROR: extract_circular_rna_features returned no features")
    avsun = _load(_AVSUN, feature_names)
    decoy = _load(_DECOY, feature_names)
    values = np.vstack((avsun, decoy))
    is_avsun = np.zeros(values.shape[0], dtype=bool)
    is_avsun[: avsun.shape[0]] = True
    _write_plot(_OUT, feature_names, values, is_avsun)
    print("wrote %s" % _OUT)
    print("avsun %d decoy %d features %d" % (avsun.shape[0], decoy.shape[0], len(feature_names)))


if __name__ == "__main__":
    main()
