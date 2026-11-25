#!/usr/bin/env python3
"""Pick control labels from big_output/controls.seqlengths.

Zetavirus contributes every label. Every other family drops lengths outside
100..6000, then keeps 100 labels spread evenly across the length-sorted list,
including the shortest and the longest. A family with fewer than 100 sequences
inside the window keeps all of them.
"""

from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
_SRC = _ROOT / "big_output" / "controls.seqlengths"
_DST = _ROOT / "big_output" / "control_labels_100_per_family.txt"
_PER_FAMILY = 100
_MIN_LEN = 100
_MAX_LEN = 6000
_ALL_FAMILIES = {"zetavirus"}


def _die(message: str) -> None:
    sys.stderr.write(message + "\n")
    sys.exit(1)


def _evenly_spaced(items: list, count: int) -> list:
    """Take count items spread across a list, including both ends."""
    if count >= len(items):
        return list(items)
    last = len(items) - 1
    chosen = [items[int(round(i * last / (count - 1)))] for i in range(count)]
    if len(chosen) != len(set(id(item) for item in chosen)):
        _die("ERROR: evenly spaced selection collided")
    return chosen


def _family(label: str) -> str:
    kind, family, rest = label.split("_", 2)
    if kind not in ("VLE", "DECOY") or not family or not rest:
        _die("ERROR: label is not VLE_family_... or DECOY_family_...: %s" % label)
    return family


def _quantile(sorted_lengths: list[int], frac: float) -> int:
    return sorted_lengths[int(round(frac * (len(sorted_lengths) - 1)))]


def _load(path: Path) -> dict[str, list[tuple[int, str]]]:
    groups: dict[str, list[tuple[int, str]]] = {}
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                _die("ERROR: %s:%d expected label and seqlength" % (path, line_number))
            label, length_text = fields
            if label in seen:
                _die("ERROR: duplicate label %s" % label)
            seen.add(label)
            try:
                length = int(length_text)
            except ValueError:
                _die("ERROR: %s:%d seqlength is not an integer" % (path, line_number))
            groups.setdefault(_family(label), []).append((length, label))
    return groups


def _summarize(
    family: str,
    total: int,
    pool: list[tuple[int, str]],
    chosen: list[tuple[int, str]],
) -> str:
    pool_lengths = [length for length, _label in pool]
    chosen_lengths = [length for length, _label in chosen]
    return "\t".join(
        [
            family,
            str(total),
            str(total - len(pool)),
            str(len(pool)),
            str(len(chosen)),
            str(chosen_lengths[0]),
            str(_quantile(chosen_lengths, 0.5)),
            str(chosen_lengths[-1]),
            str(pool_lengths[0]),
            str(_quantile(pool_lengths, 0.5)),
            str(pool_lengths[-1]),
        ]
    )


def main() -> None:
    if not _SRC.is_file():
        _die("ERROR: missing %s" % _SRC)
    groups = _load(_SRC)
    selected: list[tuple[str, int, str]] = []
    sys.stdout.write(
        "family\ttotal\toutside\tpool\tselected"
        "\tsel_min\tsel_median\tsel_max"
        "\tpool_min\tpool_median\tpool_max\n"
    )
    for family in sorted(groups):
        records = sorted(groups[family])
        if family in _ALL_FAMILIES:
            pool = records
            chosen = records
        else:
            pool = [
                record
                for record in records
                if _MIN_LEN <= record[0] <= _MAX_LEN
            ]
            if not pool:
                _die("ERROR: %s has no sequences in %d..%d" % (family, _MIN_LEN, _MAX_LEN))
            chosen = _evenly_spaced(pool, _PER_FAMILY)
        sys.stdout.write(_summarize(family, len(records), pool, chosen) + "\n")
        selected.extend((family, length, label) for length, label in chosen)

    _DST.parent.mkdir(parents=True, exist_ok=True)
    with _DST.open("w", encoding="utf-8", newline="\n") as handle:
        for _family_name, _length, label in selected:
            handle.write(label + "\n")
    sys.stdout.write("wrote %d labels to %s\n" % (len(selected), _DST))


if __name__ == "__main__":
    main()
