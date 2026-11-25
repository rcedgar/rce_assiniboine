#!/bin/bash -e
# One black PNG thumbnail per naview-only SVG. Ink is thickened in the
# bitmap. The SVG files are left unchanged.

cd "$(dirname "$0")/.."
set -o pipefail

if ! command -v convert >/dev/null; then
	echo "create_thumbnail_pngs.bash: convert is missing" >&2
	exit 1
fi

src=big_control_svgs_100_per_family_naview_only
dst=big_control_thumbnail_pngs

if [ ! -d "$src" ]; then
	echo "create_thumbnail_pngs.bash: missing $src" >&2
	exit 1
fi
mkdir -p "$dst"

export dst
find "$src" -name '*.svg' -print0 | xargs -0 -P 32 -n 1 bash -c '
set -e
svg="$1"
stem=$(basename "$svg" .svg)
out="$dst/$stem.png"
if [ -f "$out" ] && [ "$out" -nt "$svg" ]; then
	exit 0
fi
tmp="$dst/.$stem.$$.png"
trap "rm -f \"\$tmp\"" EXIT
convert -background white "$svg" \
	-alpha remove \
	-resize 480x \
	-colorspace Gray \
	-threshold 97% \
	-resize 240x \
	-threshold 50% \
	-negate -morphology Dilate Octagon:1 -negate \
	"PNG8:$tmp"
mv "$tmp" "$out"
trap - EXIT
' _

svgs=$(find "$src" -name '*.svg' | wc -l)
pngs=$(find "$dst" -name '*.png' ! -name '.*' | wc -l)
echo "wrote $pngs thumbnails from $svgs svgs to $dst"
if [ "$pngs" -ne "$svgs" ]; then
	echo "create_thumbnail_pngs.bash: expected $svgs pngs, found $pngs" >&2
	exit 1
fi
