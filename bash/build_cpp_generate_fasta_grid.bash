#!/bin/bash -e

root="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$root/build-linux"
g++ -O3 -std=c++17 -Wall -Wextra \
	-o "$root/build-linux/cpp_generate_fasta_grid" \
	"$root/cpp_generate_fasta_grid/cpp_generate_fasta_grid.cpp"
