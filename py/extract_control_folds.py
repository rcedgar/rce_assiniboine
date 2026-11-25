#!/usr/bin/env python3
"""Write selected RNAfold records as label, sequence, dot-bracket.

Input records are three lines: a FASTA header, an ACGU sequence, and a
dot-bracket line ending with the MFE in parentheses. The MFE is omitted.
Rows follow the order of the label list.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
_FOLDS = _ROOT / "big_output" / "control_RNAfold.txt"
_LABELS = _ROOT / "big_output" / "control_labels_100_per_family.txt"
_DST = _ROOT / "big_output" / "controls_folds_100_per_family.tsv"
_BASES = set("ACGTURYSWKMBDHVN")
_IUPAC = re.compile(r"^[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+$")
_STRUCTURE = re.compile(r"^([().]+)\s+\(\s*[-+]?\d+\.\d+\s*\)\s*$")


def _die(message: str) -> None:
    sys.stderr.write(message + "\n")
    sys.exit(1)


def _load_labels(path: Path) -> list[str]:
    labels: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            label = raw.strip()
            if not label:
                continue
            if label in seen:
                _die("ERROR: %s:%d duplicate label %s" % (path, line_number, label))
            seen.add(label)
            labels.append(label)
    if not labels:
        _die("ERROR: no labels in %s" % path)
    return labels


def _structure(label: str, line: str) -> str:
    match = _STRUCTURE.fullmatch(line)
    if match is None:
        _die("ERROR: %s structure line is not dot-bracket plus MFE: %s" % (label, line[:80]))
    return match.group(1)


def _extract(path: Path, wanted: set[str]) -> dict[str, tuple[str, str]]:
    found: dict[str, tuple[str, str]] = {}
    label: str | None = None
    sequence: str | None = None
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if label is None:
                if not line:
                    continue
                if not line.startswith(">"):
                    _die("ERROR: %s:%d expected a FASTA header" % (path, line_number))
                label = line[1:].split()[0]
                if not label:
                    _die("ERROR: %s:%d empty FASTA header" % (path, line_number))
                continue
            if sequence is None:
                if not line or line.startswith(">"):
                    _die("ERROR: %s:%d expected a sequence for %s" % (path, line_number, label))
                sequence = line
                continue
            if not line.endswith(")"):
                _die("ERROR: %s:%d expected a structure line for %s" % (path, line_number, label))
            if label in wanted:
                if label in found:
                    _die("ERROR: duplicate RNAfold record for %s" % label)
                if _IUPAC.fullmatch(sequence) is None:
                    odd = "".join(sorted(set(ch for ch in sequence if ch.upper() not in _BASES)))
                    _die("ERROR: %s sequence has unexpected bases %s" % (label, odd))
                dot_bracket = _structure(label, line)
                if len(dot_bracket) != len(sequence):
                    _die(
                        "ERROR: %s sequence length %d differs from dot-bracket length %d"
                        % (label, len(sequence), len(dot_bracket))
                    )
                found[label] = (sequence, dot_bracket)
                if len(found) == len(wanted):
                    return found
            label = None
            sequence = None
    if label is not None:
        _die("ERROR: %s ended inside the record for %s" % (path, label))
    return found


def main() -> None:
    if not _FOLDS.is_file():
        _die("ERROR: missing %s" % _FOLDS)
    if not _LABELS.is_file():
        _die("ERROR: missing %s" % _LABELS)
    labels = _load_labels(_LABELS)
    found = _extract(_FOLDS, set(labels))
    missing = [label for label in labels if label not in found]
    if missing:
        preview = ", ".join(missing[:8])
        _die("ERROR: %d labels missing from %s: %s" % (len(missing), _FOLDS, preview))

    _DST.parent.mkdir(parents=True, exist_ok=True)
    lengths: list[int] = []
    ambiguous = 0
    lowercase = 0
    with _DST.open("w", encoding="utf-8", newline="\n") as handle:
        for label in labels:
            sequence, dot_bracket = found[label]
            if "\t" in label or "\t" in sequence or "\t" in dot_bracket:
                _die("ERROR: tab in a field for %s" % label)
            handle.write("%s\t%s\t%s\n" % (label, sequence, dot_bracket))
            lengths.append(len(sequence))
            if any(base not in "ACGU" for base in sequence):
                ambiguous += 1
            if any(base.islower() for base in sequence):
                lowercase += 1
    lengths.sort()
    sys.stdout.write(
        "wrote %d rows to %s; L min %d median %d max %d; "
        "outside ACGU: %d; containing lowercase: %d\n"
        % (
            len(labels),
            _DST,
            lengths[0],
            lengths[len(lengths) // 2],
            lengths[-1],
            ambiguous,
            lowercase,
        )
    )


if __name__ == "__main__":
    main()
