#!/bin/bash -e

root="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$root/big_fasta_grid"
"$root/build-linux/cpp_generate_fasta_grid" \
	"$root/output/lcomp_bin_ranges.tsv" \
	"$root/big_fasta_grid" \
	100
