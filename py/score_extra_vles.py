#!/usr/bin/env python3
"""Score 50 more VLE sequences and append them to the control Z table.

Candidates are VLE records in big_output/controls.fa that are ACGT, inside the
RNAfold grid length range, and not already present in controls_z_scores.tsv.
The 50 sequences are spread evenly across that list, including both ends.
"""

from __future__ import annotations

import sys
from pathlib import Path

import correlate_measured_predicted_Z as corr
import mfe_grid_interp

_ROOT = Path(__file__).resolve().parent.parent
_WANTED = 50


def _read_labels(path: Path) -> set[str]:
    labels: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        index = header.index("label")
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) > index and fields[index]:
                labels.add(fields[index])
    return labels


def _vle_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    seq_id = ""
    parts: list[str] = []

    def finish() -> None:
        if seq_id.startswith("VLE_") and parts:
            records.append((seq_id, "".join(parts)))

    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                seq_id = header.split()[0] if header else ""
                parts = []
            elif seq_id.startswith("VLE_"):
                parts.append(line.replace(" ", "").replace("\t", ""))
    finish()
    return records


def main() -> None:
    outdir = _ROOT / "correlate_measured_predicted_Z"
    controls_path = outdir / "controls_z_scores.tsv"
    extra_path = outdir / "vle_extra50_z_scores.tsv"
    already = _read_labels(controls_path)
    table_rows, lengths, k = mfe_grid_interp.load_table(str(_ROOT / "output" / "rnafold_grid_mean_var.tsv"))
    candidates = []
    for seq_id, seq in _vle_records(_ROOT / "big_output" / "controls.fa"):
        if seq_id in already:
            continue
        if corr._sequence_can_score(seq, lengths):
            candidates.append((seq_id, seq))
    if len(candidates) < _WANTED:
        corr._die(f"ERROR: only {len(candidates)} new scorable VLEs")
    chosen = corr._evenly_spaced(candidates, _WANTED)
    lengths_chosen = [len(seq) for _seq_id, seq in chosen]
    print(
        "scoring %d of %d new scorable VLEs; L min %d median %d max %d"
        % (
            len(chosen),
            len(candidates),
            min(lengths_chosen),
            lengths_chosen[len(lengths_chosen) // 2],
            max(lengths_chosen),
        ),
        file=sys.stderr,
        flush=True,
    )
    cpu_threads = 64
    rows = []
    with extra_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(corr._COLUMNS) + "\n")
        for index, (seq_id, seq) in enumerate(chosen, start=1):
            row, detail = corr.score_sequence(
                seq_id,
                seq,
                seq_id,
                table_rows,
                lengths,
                k,
                100,
                1,
                cpu_threads,
            )
            handle.write(corr._format_row(row) + "\n")
            handle.flush()
            rows.append(row)
            corr._progress(index, len(chosen), row, detail)
    ok = [row for row in rows if row["status"] == "ok"]
    if len(ok) < _WANTED:
        corr._die(f"ERROR: only {len(ok)} of {_WANTED} VLEs scored")
    with controls_path.open("a", encoding="utf-8", newline="\n") as handle:
        for row in ok:
            handle.write(corr._format_row(row) + "\n")
    print("appended %d VLEs to %s" % (len(ok), controls_path), file=sys.stderr)


if __name__ == "__main__":
    main()
