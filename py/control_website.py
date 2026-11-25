#!/usr/bin/env python3
"""Build a flat HTML site for the control-sequence figures.

index.html shows one thumbnail per family. A family is named
VLE_<family> or DECOY_<family>. VLE families come first, then DECOY
families start on a new row. Each card is the same width. VLE cards
have a 90% transparent green background and DECOY cards a 90%
transparent red background. Each family page shows every
sequence in that family, longest first, as a 5-column grid that fills
the viewport. Thumbnails are PNG files. Each sequence page wraps the
full SVG. Prev, Next, Up, and Home sit left of the title. A button
with nowhere to go stays put and is gray. Arrow keys and Home follow
those buttons.

The index thumbnail for a family is chosen with random.Random(1).
"""

from __future__ import annotations

import html
import os
import random
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
_LABELS = _ROOT / "big_output" / "control_labels_100_per_family.txt"
_LENGTHS = _ROOT / "big_output" / "controls.seqlengths"
_FULL = _ROOT / "big_control_svgs_100_per_family"
_THUMB = _ROOT / "big_control_thumbnail_pngs"
_HTML = _ROOT / "big_control_website"
_PER_FAMILY = 100
_COLUMNS = 5
_FAMILY_CHARS = set(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._"
)
_KINDS = ("VLE", "DECOY")

STYLE = """
body { margin: 0; font-family: sans-serif; }
h1 { font-size: 1.25rem; font-weight: 600; margin: 0; min-width: 0; }
.head {
  display: flex;
  align-items: center;
  gap: 12px;
  margin: 16px 16px 8px;
}
nav.jump { display: flex; flex: 0 0 auto; gap: 6px; }
nav.jump a,
nav.jump span {
  font-size: 0.8rem;
  line-height: 1.2;
  padding: 2px 8px;
  border: 1px solid #888;
  border-radius: 3px;
  background: #f4f4f4;
  color: #0645ad;
  text-decoration: none;
}
nav.jump a:hover { background: #e4e4e4; }
nav.jump span {
  color: #888;
  border-color: #ccc;
  background: #eee;
  cursor: default;
}
p.crumb { margin: 0 16px 8px; }
p.crumb a { color: #0645ad; text-decoration: none; }
p.crumb a:hover { text-decoration: underline; }
.families {
  display: flex;
  flex-wrap: wrap;
  gap: 12px 16px;
  margin: 0 16px 24px;
  align-items: end;
}
.families .rowbreak {
  flex-basis: 100%;
  height: 0;
}
.families a {
  display: block;
  width: 160px;
  border: 1px solid #888;
  padding: 3px;
  color: inherit;
  text-decoration: none;
}
.families a.vle,
.thumbs.vle a { background: rgba(0, 128, 0, 0.1); }
.families a.decoy,
.thumbs.decoy a { background: rgba(204, 0, 0, 0.1); }
.families img {
  width: 100%;
  height: auto;
  display: block;
  mix-blend-mode: multiply;
}
.families .name {
  display: block;
  margin-top: 2px;
  font-size: 13px;
  line-height: 1.2;
  text-align: center;
  color: #222;
}
body.gridpage {
  height: 100vh;
  display: flex;
  flex-direction: column;
  overflow: hidden;
}
body.gridpage .head { margin: 8px 16px 4px; }
body.gridpage .head h1 { font-size: 1rem; }
.thumbs {
  flex: 1 1 auto;
  min-height: 0;
  display: grid;
  grid-template-columns: repeat(5, minmax(0, 1fr));
  gap: 4px;
  margin: 0 8px 8px;
}
.thumbs a {
  min-height: 0;
  min-width: 0;
  display: flex;
  align-items: center;
  justify-content: center;
  overflow: hidden;
  border: 1px solid #888;
  padding: 2px;
}
.thumbs img {
  max-width: 100%;
  max-height: 100%;
  object-fit: contain;
  mix-blend-mode: multiply;
}
body.svg {
  height: 100vh;
  display: flex;
  flex-direction: column;
  overflow: hidden;
}
body.svg .head { margin: 8px 16px 4px; }
body.svg .head h1 {
  font-size: 1rem;
  overflow-wrap: anywhere;
}
body.svg p.crumb { margin: 0 16px 4px; }
body.svg .view {
  flex: 1 1 auto;
  min-height: 0;
  display: flex;
  align-items: flex-start;
  justify-content: flex-start;
  overflow: auto;
}
body.svg .view img {
  display: block;
  max-width: 100%;
  max-height: 100%;
  width: auto;
  height: auto;
  object-fit: contain;
}
""".strip()

NAV_SCRIPT = """
(function () {
  function go(id) {
    var el = document.getElementById(id);
    if (el && el.tagName === "A" && el.getAttribute("href")) {
      location.href = el.href;
    }
  }
  document.addEventListener("keydown", function (event) {
    if (event.defaultPrevented) return;
    if (event.altKey || event.ctrlKey || event.metaKey) return;
    var target = event.target;
    var tag = (target && target.tagName) || "";
    if (tag === "INPUT" || tag === "TEXTAREA" || (target && target.isContentEditable)) {
      return;
    }
    if (event.key === "ArrowLeft") {
      event.preventDefault();
      go("prev");
    } else if (event.key === "ArrowRight") {
      event.preventDefault();
      go("next");
    } else if (event.key === "ArrowUp") {
      event.preventDefault();
      go("up");
    } else if (event.key === "ArrowDown") {
      event.preventDefault();
      go("down");
    } else if (event.key === "Home") {
      event.preventDefault();
      go("home");
    }
  });
})();
""".strip()


class Seq:
    def __init__(self, label: str, family: str, length: int, stem: str):
        self.label = label
        self.family = family
        self.length = length
        self.stem = stem


def _die(message: str) -> None:
    sys.stderr.write(message + "\n")
    sys.exit(1)


def _is_windows_device_name(stem: str) -> bool:
    folded = stem.upper()
    if folded in ("CON", "PRN", "AUX", "NUL"):
        return True
    if len(folded) == 4 and folded[:3] in ("COM", "LPT") and folded[3] in "123456789":
        return True
    return False


def _stem_from_label(label: str, line_no: int) -> str:
    """Match rce_naview stem_from_label: keep [A-Za-z0-9._], cut at 80."""
    chars = []
    for char in label:
        keep = (
            ("A" <= char <= "Z")
            or ("a" <= char <= "z")
            or ("0" <= char <= "9")
            or char in "._"
        )
        chars.append(char if keep else "_")
    stem = "".join(chars)
    if len(stem) > 80:
        stem = stem[:80]
    if stem in ("", ".", "..") or _is_windows_device_name(stem):
        stem = "record_%d" % line_no
    return stem


def _family_of(label: str) -> str:
    parts = label.split("_", 2)
    if len(parts) != 3 or parts[0] not in _KINDS or not parts[1] or not parts[2]:
        _die("ERROR: label is not VLE_family_... or DECOY_family_...: %s" % label)
    family = parts[1]
    if any(char not in _FAMILY_CHARS for char in family):
        _die("ERROR: family name is not a safe file stem: %s" % family)
    return "%s_%s" % (parts[0], family)


def _family_sort_key(family: str) -> tuple[int, str]:
    kind, rest = family.split("_", 1)
    return (0 if kind == "VLE" else 1, rest)


def _load_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                _die("ERROR: %s:%d expected label and seqlength" % (path, line_number))
            label, length_text = fields
            if label in lengths:
                _die("ERROR: duplicate label %s" % label)
            try:
                length = int(length_text)
            except ValueError:
                _die("ERROR: %s:%d seqlength is not an integer" % (path, line_number))
            lengths[label] = length
    return lengths


def _load_sequences(labels_path: Path, lengths: dict[str, int]) -> list[Seq]:
    sequences: list[Seq] = []
    seen: set[str] = set()
    stems: dict[str, str] = {}
    with labels_path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            label = raw.strip()
            if not label:
                continue
            if label in seen:
                _die("ERROR: %s:%d duplicate label %s" % (labels_path, line_number, label))
            seen.add(label)
            if label not in lengths:
                _die("ERROR: no seqlength for %s" % label)
            family = _family_of(label)
            stem = _stem_from_label(label, line_number)
            previous = stems.get(stem)
            if previous is not None:
                _die("ERROR: %s and %s both write %s.svg" % (previous, label, stem))
            stems[stem] = label
            full = _FULL / ("%s.svg" % stem)
            thumb = _THUMB / ("%s.png" % stem)
            if not full.is_file():
                _die("ERROR: missing %s" % full)
            if not thumb.is_file():
                _die("ERROR: missing %s" % thumb)
            sequences.append(Seq(label, family, lengths[label], stem))
    if not sequences:
        _die("ERROR: no labels in %s" % labels_path)
    return sequences


def _group(sequences: list[Seq]) -> dict[str, list[Seq]]:
    groups: dict[str, list[Seq]] = {}
    for seq in sequences:
        groups.setdefault(seq.family, []).append(seq)
    for family, members in groups.items():
        if len(members) > _PER_FAMILY:
            _die("ERROR: %s has %d sequences" % (family, len(members)))
        members.sort(key=lambda item: (-item.length, item.label))
    return groups


def _family_page(family: str) -> str:
    return "%s_index.html" % family


def _sequence_page(stem: str) -> str:
    return "%s.html" % stem


def _document(title: str, body: str, body_class: str | None = None) -> str:
    body_attr = ' class="%s"' % body_class if body_class else ""
    return (
        "<!DOCTYPE html>\n"
        "<html>\n"
        "<head>\n"
        '<meta charset="utf-8">\n'
        "<title>%s</title>\n"
        "<style>\n%s\n</style>\n"
        "</head>\n"
        "<body%s>\n"
        "%s\n"
        "</body>\n"
        "</html>\n"
    ) % (html.escape(title), STYLE, body_attr, body)


def _write_html(htmldir: Path, name: str, text: str) -> None:
    path = htmldir / name
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _clear_html(htmldir: Path) -> None:
    if not htmldir.is_dir():
        return
    for path in htmldir.glob("*.html"):
        path.unlink()


def _svg_src(svg_path: Path, htmldir: Path) -> str:
    stamp = int(svg_path.stat().st_mtime)
    relative = Path(os.path.relpath(svg_path, htmldir)).as_posix()
    return "%s?v=%d" % (relative, stamp)


def _nav_button(element_id: str, href: str | None, label: str, shortcut_name: str) -> str:
    tip = html.escape(shortcut_name, quote=True)
    text = html.escape(label)
    if not href:
        return '<span title="%s" aria-disabled="true">%s</span>' % (tip, text)
    return '<a id="%s" href="%s" title="%s">%s</a>' % (
        html.escape(element_id, quote=True),
        html.escape(href, quote=True),
        tip,
        text,
    )


def _shortcut(element_id: str, href: str | None) -> str:
    if not href:
        return ""
    return '<a id="%s" href="%s" hidden></a>' % (
        html.escape(element_id, quote=True),
        html.escape(href, quote=True),
    )


def _neighbors(names: list[str], index: int) -> tuple[str | None, str | None]:
    prev_href = names[index - 1] if index else None
    next_href = names[index + 1] if index + 1 < len(names) else None
    return prev_href, next_href


def _page_nav(
    up_href: str | None,
    prev_href: str | None,
    next_href: str | None,
    home_href: str | None = "index.html",
    down_href: str | None = None,
) -> str:
    buttons = [
        _nav_button("prev", prev_href, "Prev", "Left arrow"),
        _nav_button("next", next_href, "Next", "Right arrow"),
        _nav_button("up", up_href, "Up", "Up arrow"),
        _nav_button("home", home_href, "Home", "Home"),
    ]
    parts = ['<nav class="jump">%s</nav>' % "".join(buttons)]
    down = _shortcut("down", down_href)
    if down:
        parts.append(down)
    parts.append("<script>\n%s\n</script>" % NAV_SCRIPT)
    return "\n".join(parts)


def _rows(count: int) -> int:
    return (count + _COLUMNS - 1) // _COLUMNS


def _build(htmldir: Path, groups: dict[str, list[Seq]]) -> None:
    htmldir.mkdir(parents=True, exist_ok=True)
    _clear_html(htmldir)
    families = sorted(groups, key=_family_sort_key)
    family_hrefs = [_family_page(family) for family in families]
    rng = random.Random(1)
    picks = {family: rng.choice(groups[family]) for family in families}

    cards = []
    started_decoy = False
    for family in families:
        if family.startswith("DECOY_") and not started_decoy:
            cards.append('<div class="rowbreak"></div>')
            started_decoy = True
        kind = "vle" if family.startswith("VLE_") else "decoy"
        pick = picks[family]
        src = _svg_src(_THUMB / ("%s.png" % pick.stem), htmldir)
        cards.append(
            '<a class="%s" href="%s"><img src="%s" alt="%s"><span class="name">%s</span></a>'
            % (
                kind,
                html.escape(_family_page(family), quote=True),
                html.escape(src, quote=True),
                html.escape(family),
                html.escape(family),
            )
        )
    index_body = "\n".join(
        [
            '<div class="head">',
            _page_nav(None, None, None, home_href=None, down_href=family_hrefs[0]),
            "<h1>Controls</h1>",
            "</div>",
            '<div class="families">%s</div>' % "".join(cards),
        ]
    )
    _write_html(htmldir, "index.html", _document("Controls", index_body))

    for family_index, family in enumerate(families):
        members = groups[family]
        prev_href, next_href = _neighbors(family_hrefs, family_index)
        cells = []
        for seq in members:
            src = _svg_src(_THUMB / ("%s.png" % seq.stem), htmldir)
            cells.append(
                '<a href="%s" title="%s"><img src="%s" alt="%s"></a>'
                % (
                    html.escape(_sequence_page(seq.stem), quote=True),
                    html.escape(seq.label, quote=True),
                    html.escape(src, quote=True),
                    html.escape(seq.label),
                )
            )
        family_body = "\n".join(
            [
                '<div class="head">',
                _page_nav(
                    "index.html",
                    prev_href,
                    next_href,
                    down_href=_sequence_page(members[0].stem),
                ),
                "<h1>%s</h1>" % html.escape(family),
                "</div>",
                '<div class="thumbs %s" style="grid-template-rows: repeat(%d, minmax(0, 1fr))">%s</div>'
                % (
                    "vle" if family.startswith("VLE_") else "decoy",
                    _rows(len(members)),
                    "".join(cells),
                ),
            ]
        )
        _write_html(
            htmldir,
            _family_page(family),
            _document(family, family_body, body_class="gridpage"),
        )

        sequence_hrefs = [_sequence_page(seq.stem) for seq in members]
        up_href = _family_page(family)
        for seq_index, seq in enumerate(members):
            seq_prev, seq_next = _neighbors(sequence_hrefs, seq_index)
            src = _svg_src(_FULL / ("%s.svg" % seq.stem), htmldir)
            crumb = '<p class="crumb"><a href="%s">%s</a></p>' % (
                html.escape(up_href, quote=True),
                html.escape(family),
            )
            seq_body = "\n".join(
                [
                    '<div class="head">',
                    _page_nav(up_href, seq_prev, seq_next),
                    "<h1>%s</h1>" % html.escape(seq.label),
                    "</div>",
                    crumb,
                    '<div class="view">',
                    '<img src="%s" alt="%s">'
                    % (html.escape(src, quote=True), html.escape(seq.label)),
                    "</div>",
                ]
            )
            _write_html(
                htmldir,
                _sequence_page(seq.stem),
                _document(seq.label, seq_body, body_class="svg"),
            )


def main() -> None:
    if not _LABELS.is_file():
        _die("ERROR: missing %s" % _LABELS)
    if not _LENGTHS.is_file():
        _die("ERROR: missing %s" % _LENGTHS)
    lengths = _load_lengths(_LENGTHS)
    sequences = _load_sequences(_LABELS, lengths)
    groups = _group(sequences)
    _build(_HTML, groups)
    sys.stdout.write(
        "wrote %d sequence pages in %d families to %s\n"
        % (len(sequences), len(groups), _HTML)
    )


if __name__ == "__main__":
    main()
