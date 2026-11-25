#!/bin/bash -e

root="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$root/build-linux"
g++ -O3 -std=c++17 -Wall -Wextra -pthread \
	-o "$root/build-linux/rce_naview" \
	"$root/rce_naview/main.cpp" \
	"$root/rce_naview/naview.cpp" \
	"$root/rce_naview/svg_page.cpp"
