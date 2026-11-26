#!/bin/bash -e
# Composition-colored control site. Each base is painted from a 20-base
# window. Rainbow SVGs and the default site are left in place.

cd "$(dirname "$0")/.."
set -o pipefail

window=20
folds=big_output/controls_folds_100_per_family.tsv
full=big_control_svgs_100_per_family_w${window}
nav=big_control_svgs_100_per_family_w${window}_naview_only
png=big_control_thumbnail_pngs_w${window}
naview=build-linux/rce_naview

if [ ! -x "$naview" ]; then
	echo "control_website_w20.bash: missing $naview" >&2
	exit 1
fi
if [ ! -f "$folds" ]; then
	echo "control_website_w20.bash: missing $folds" >&2
	exit 1
fi
if ! command -v convert >/dev/null; then
	echo "control_website_w20.bash: convert is missing" >&2
	exit 1
fi

mkdir -p "$full" "$nav" "$png"

"$naview" \
	--circular \
	--color-window "$window" \
	--threads 32 \
	--out-dir "$full" \
	"$folds"

"$naview" \
	--circular \
	--naview-only \
	--color-window "$window" \
	--threads 32 \
	--out-dir "$nav" \
	"$folds"

export png
find "$nav" -name '*.svg' -print0 | xargs -0 -P 32 -n 1 bash -c '
set -e
svg="$1"
stem=$(basename "$svg" .svg)
out="$png/$stem.png"
if [ -f "$out" ] && [ "$out" -nt "$svg" ]; then
	exit 0
fi
tmp="$png/.$stem.$$.png"
trap "rm -f \"\$tmp\"" EXIT
convert -background white "$svg" \
	-alpha remove \
	-resize 240x \
	"PNG32:$tmp"
mv "$tmp" "$out"
trap - EXIT
' _

svgs=$(find "$nav" -name '*.svg' | wc -l)
pngs=$(find "$png" -name '*.png' ! -name '.*' | wc -l)
echo "wrote $pngs thumbnails from $svgs svgs to $png"
if [ "$pngs" -ne "$svgs" ]; then
	echo "control_website_w20.bash: expected $svgs pngs, found $pngs" >&2
	exit 1
fi

python3 py/control_website.py --color-window "$window"
