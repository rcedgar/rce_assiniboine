#!/bin/bash -e
# VLE_avsunviroid thumbnails and length-matched non-avsun decoys.
# Sequence length minL..maxL is taken from the avsun rows.
# Each other VLE_<family> or DECOY_<family> keeps at most 1000 rows,
# evenly spaced by sequence length, including the shortest and longest.
# Avsun rows are dropped when pct_bulge > 10, pH < 17, or
# bulge_to_stem_ratio > 0.14.

cd "$(dirname "$0")/.."
set -o pipefail

bash bash/build_rce_naview.bash

src=big_output/controls_label_seq_ss.tsv
out=big_avsum_thumbs
prefix=VLE_avsunviroid

if [ ! -f "$src" ]; then
	echo "avsun_thumbnails.bash: missing $src" >&2
	exit 1
fi

rm -rf "$out/decoy_thumbnails" "$out/misc/decoy_parts"
mkdir -p "$out/misc" "$out/avsun_thumbnails" "$out/decoy_thumbnails" "$out/misc/decoy_parts"

avsun="$out/misc/avsun.tsv"
decoy="$out/misc/decoy.tsv"
range="$out/misc/length_range.txt"
: > "$avsun"
: > "$decoy"

read -r minL maxL avsun_n < <(awk -F '\t' -v pre="$prefix" -v out="$avsun" '
BEGIN { minL = 0; maxL = 0; n = 0 }
index($1, pre) == 1 {
	L = length($2)
	if (n == 0 || L < minL) minL = L
	if (n == 0 || L > maxL) maxL = L
	n++
	print > out
	fflush(out)
}
END { print minL, maxL, n }
' "$src")

if [ "$avsun_n" -lt 1 ]; then
	echo "avsun_thumbnails.bash: no $prefix rows" >&2
	exit 1
fi

read -r avsun_kept avsun_discarded < <(
	/home/bob/miniconda3/envs/null_model/bin/python -u py/avsun_feature_distributions.py --filter-only
)
avsun_n=$avsun_kept
echo "kept $avsun_kept avsun, discarded $avsun_discarded"

export AVSUN_PREFIX="$prefix" AVSUN_MIN="$minL" AVSUN_MAX="$maxL"
export AVSUN_SRC="$src" AVSUN_DECOY="$decoy" AVSUN_CAP=1000
read -r decoy_pool decoy_n decoy_families < <(python3 - <<'PY'
import os
import sys

prefix = os.environ["AVSUN_PREFIX"]
min_l = int(os.environ["AVSUN_MIN"])
max_l = int(os.environ["AVSUN_MAX"])
cap = int(os.environ["AVSUN_CAP"])
src = os.environ["AVSUN_SRC"]
dst = os.environ["AVSUN_DECOY"]


def spaced_indices(n, k):
    if k >= n:
        return list(range(n))
    if k == 1:
        return [0]
    idxs = []
    prev = -1
    for i in range(k):
        idx = (i * (n - 1) + (k - 1) // 2) // (k - 1)
        if idx <= prev:
            idx = prev + 1
        remain = k - 1 - i
        max_idx = n - 1 - remain
        if idx > max_idx:
            idx = max_idx
        idxs.append(idx)
        prev = idx
    return idxs


def family_of(label):
    parts = label.split("_", 2)
    if len(parts) != 3 or parts[0] not in ("VLE", "DECOY") or not parts[1] or not parts[2]:
        sys.exit("avsun_thumbnails.bash: label is not VLE_family_... or DECOY_family_...: " + label)
    return parts[0] + "_" + parts[1]


rows = []
with open(src, "r", encoding="utf-8", newline="") as handle:
    for raw in handle:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            sys.exit("avsun_thumbnails.bash: expected 3 columns")
        label, seq, ss = fields[0], fields[1], fields[2]
        if label.startswith(prefix):
            continue
        if len(seq) != len(ss):
            sys.exit("avsun_thumbnails.bash: sequence and structure lengths differ: " + label)
        length = len(seq)
        if length < min_l or length > max_l:
            continue
        rows.append((family_of(label), length, line))

groups = {}
for index, row in enumerate(rows):
    groups.setdefault(row[0], []).append(index)

keep = []
for indexes in groups.values():
    ordered = sorted(indexes, key=lambda index: (rows[index][1], index))
    for pick in spaced_indices(len(ordered), min(cap, len(ordered))):
        keep.append(ordered[pick])
keep.sort()
if len(keep) != len(set(keep)):
    sys.exit("avsun_thumbnails.bash: duplicate decoy picks")

with open(dst, "w", encoding="utf-8", newline="\n") as handle:
    for index in keep:
        handle.write(rows[index][2] + "\n")

print(len(rows), len(keep), len(groups))
PY
)

printf 'minL %s\nmaxL %s\navsun %s\navsun_discarded %s\ndecoy_pool %s\ndecoy %s\ncap %s\nfamilies %s\n' \
	"$minL" "$maxL" "$avsun_n" "$avsun_discarded" "$decoy_pool" "$decoy_n" "1000" "$decoy_families" > "$range"
echo "minL $minL maxL $maxL avsun $avsun_n decoy_pool $decoy_pool decoy $decoy_n families $decoy_families"

nav=build-linux/rce_naview
avsun_lines=$(wc -l < "$avsun")
avsun_png_now=$(find "$out/avsun_thumbnails" -name '*.png' ! -name '.*' | wc -l)
if [ "$avsun_png_now" -eq "$avsun_lines" ]; then
	echo "avsun thumbnails already complete ($avsun_png_now)"
else
	rm -rf "$out/avsun_thumbnails"
	mkdir -p "$out/avsun_thumbnails"
	echo "drawing avsun"
	"$nav" --circular --png 200,100 --threads 32 \
		--out-dir "$out/avsun_thumbnails" "$avsun"
fi

split -d -a 3 -l 4000 --additional-suffix=.tsv \
	"$decoy" "$out/misc/decoy_parts/part_"

for part in "$out/misc/decoy_parts"/part_*.tsv; do
	echo "drawing $part"
	"$nav" --circular --png 200,100 --threads 32 \
		--out-dir "$out/decoy_thumbnails" "$part"
done

avsun_lines=$(wc -l < "$avsun")
decoy_lines=$(wc -l < "$decoy")
avsun_png=$(find "$out/avsun_thumbnails" -name '*.png' ! -name '.*' | wc -l)
decoy_png=$(find "$out/decoy_thumbnails" -name '*.png' ! -name '.*' | wc -l)

echo "wrote $avsun_png avsun and $decoy_png decoy thumbnails"
if [ "$avsun_png" -ne "$avsun_lines" ] || [ "$decoy_png" -ne "$decoy_lines" ]; then
	echo "avsun_thumbnails.bash: expected $avsun_lines and $decoy_lines pngs, found $avsun_png and $decoy_png" >&2
	exit 1
fi
