#!/bin/bash -e

root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
bash "$root/bash/build_rce_naview.bash"
g++ -O2 -std=c++17 -Wall -Wextra \
	-o "$root/build-linux/rce_naview_test" \
	"$root/rce_naview/naview_test.cpp" \
	"$root/rce_naview/naview.cpp" \
	"$root/rce_naview/svg_page.cpp" \
	"$root/rce_naview/png_page.cpp"
"$root/build-linux/rce_naview_test"

mkdir -p "$root/build-linux/naview_test_out"
printf 'hairpin\tGGGAAACCC\t(((...)))\n' > "$root/build-linux/naview_test_out/in.tsv"
"$root/build-linux/rce_naview" --out-dir "$root/build-linux/naview_test_out" \
	"$root/build-linux/naview_test_out/in.tsv"
test -s "$root/build-linux/naview_test_out/hairpin.svg"
grep -q '#c62828' "$root/build-linux/naview_test_out/hairpin.svg"
"$root/build-linux/rce_naview" --rainbow --out-dir "$root/build-linux/naview_test_out" \
	"$root/build-linux/naview_test_out/in.tsv"
grep -q 'hsl(0.00, 72%, 40%)' "$root/build-linux/naview_test_out/hairpin.svg"
"$root/build-linux/rce_naview" --colors 'lime, rgb(1, 2, 3)' \
	--out-dir "$root/build-linux/naview_test_out" \
	"$root/build-linux/naview_test_out/in.tsv"
grep -q 'lime' "$root/build-linux/naview_test_out/hairpin.svg"
grep -q 'rgb(1, 2, 3)' "$root/build-linux/naview_test_out/hairpin.svg"
if "$root/build-linux/rce_naview" --rainbow --colors red \
	--out-dir "$root/build-linux/naview_test_out" \
	"$root/build-linux/naview_test_out/in.tsv"; then
	echo "expected --rainbow and --colors together to fail"
	exit 1
fi
printf 'guanine\tGGGGGGGGG\t.........\n' > "$root/build-linux/naview_test_out/guanine.tsv"
"$root/build-linux/rce_naview" --color-window 20 --naview-only \
	--out-dir "$root/build-linux/naview_test_out" \
	"$root/build-linux/naview_test_out/guanine.tsv"
grep -q '#ffa360' "$root/build-linux/naview_test_out/guanine.svg"
if grep -q 'hsl(' "$root/build-linux/naview_test_out/guanine.svg"; then
	echo "expected --color-window to replace rainbow colors"
	exit 1
fi
if "$root/build-linux/rce_naview" --color-window 20 --rainbow \
	--out-dir "$root/build-linux/naview_test_out" \
	"$root/build-linux/naview_test_out/guanine.tsv"; then
	echo "expected --color-window and --rainbow together to fail"
	exit 1
fi
if "$root/build-linux/rce_naview" --colors 'notacolor' \
	--out-dir "$root/build-linux/naview_test_out" \
	"$root/build-linux/naview_test_out/in.tsv"; then
	echo "expected an invalid color to fail"
	exit 1
fi

printf 'stemA\tGGGAAACCC\t(((...)))\nstemB\tAAAAAAA\t((...))\n' \
	> "$root/build-linux/naview_test_out/threads.tsv"
mkdir -p "$root/build-linux/naview_test_out/thr1" "$root/build-linux/naview_test_out/thr2"
"$root/build-linux/rce_naview" --threads 1 \
	--out-dir "$root/build-linux/naview_test_out/thr1" \
	"$root/build-linux/naview_test_out/threads.tsv"
"$root/build-linux/rce_naview" --threads 2 \
	--out-dir "$root/build-linux/naview_test_out/thr2" \
	"$root/build-linux/naview_test_out/threads.tsv"
cmp "$root/build-linux/naview_test_out/thr1/stemA.svg" \
	"$root/build-linux/naview_test_out/thr2/stemA.svg"
cmp "$root/build-linux/naview_test_out/thr1/stemB.svg" \
	"$root/build-linux/naview_test_out/thr2/stemB.svg"
printf 'ok\tAA\t..\nbad\tA\t(\n' > "$root/build-linux/naview_test_out/threads_bad.tsv"
rm -f "$root/build-linux/naview_test_out/thr2/ok.svg"
if "$root/build-linux/rce_naview" --threads 2 \
	--out-dir "$root/build-linux/naview_test_out/thr2" \
	"$root/build-linux/naview_test_out/threads_bad.tsv"; then
	echo "expected a bad structure to fail under --threads 2"
	exit 1
fi
test ! -e "$root/build-linux/naview_test_out/thr2/ok.svg"
if "$root/build-linux/rce_naview" --threads 0 \
	--out-dir "$root/build-linux/naview_test_out/thr1" \
	"$root/build-linux/naview_test_out/threads.tsv"; then
	echo "expected --threads 0 to fail"
	exit 1
fi

printf 'bad\tA\t(\n' > "$root/build-linux/naview_test_out/bad.tsv"
if "$root/build-linux/rce_naview" --out-dir "$root/build-linux/naview_test_out" \
	"$root/build-linux/naview_test_out/bad.tsv"; then
	echo "expected a bad structure to fail"
	exit 1
fi

mkdir -p "$root/build-linux/naview_test_out/png"
rm -f "$root/build-linux/naview_test_out/png/hairpin.svg" \
	"$root/build-linux/naview_test_out/png/hairpin.png"
"$root/build-linux/rce_naview" --png 48,32 \
	--out-dir "$root/build-linux/naview_test_out/png" \
	"$root/build-linux/naview_test_out/in.tsv"
test -s "$root/build-linux/naview_test_out/png/hairpin.png"
test ! -e "$root/build-linux/naview_test_out/png/hairpin.svg"
python3 - "$root/build-linux/naview_test_out/png/hairpin.png" <<'PY'
import struct
import sys
import zlib

data = open(sys.argv[1], "rb").read()
if data[:8] != b"\x89PNG\r\n\x1a\n":
    sys.exit("bad signature")
pos = 8
idat = b""
width = height = depth = color = None
while pos + 12 <= len(data):
    length = struct.unpack(">I", data[pos:pos + 4])[0]
    kind = data[pos + 4:pos + 8]
    pos += 8
    chunk = data[pos:pos + length]
    pos += length + 4
    if kind == b"IHDR":
        width, height, depth, color = struct.unpack(">IIBB", chunk[:10])
    elif kind == b"IDAT":
        idat += chunk
    elif kind == b"IEND":
        break
if (width, height, depth, color) != (48, 32, 1, 0):
    sys.exit("bad ihdr")
raw = zlib.decompress(idat)
stride = (48 + 7) // 8
row_bytes = 1 + stride
if len(raw) != 32 * row_bytes:
    sys.exit("bad image length")
black = False
white = False
for y in range(32):
    row = raw[y * row_bytes:(y + 1) * row_bytes]
    if row[0] != 0:
        sys.exit("bad filter")
    packed = row[1:]
    for x in range(48):
        bit = (packed[x // 8] >> (7 - (x % 8))) & 1
        if bit == 0:
            black = True
        else:
            white = True
if not black or not white:
    sys.exit("expected black and white pixels")
if ((raw[1] >> 7) & 1) != 1:
    sys.exit("expected a white corner")
PY
for bad in 0,10 10 10,10,10 no 8193,1; do
	if "$root/build-linux/rce_naview" --png "$bad" \
		--out-dir "$root/build-linux/naview_test_out/png" \
		"$root/build-linux/naview_test_out/in.tsv"; then
		echo "expected --png $bad to fail"
		exit 1
	fi
done
