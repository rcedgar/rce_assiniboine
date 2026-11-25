#!/bin/bash -e

root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
bash "$root/bash/build_rce_naview.bash"
g++ -O2 -std=c++17 -Wall -Wextra \
	-o "$root/build-linux/rce_naview_test" \
	"$root/rce_naview/naview_test.cpp" \
	"$root/rce_naview/naview.cpp" \
	"$root/rce_naview/svg_page.cpp"
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
