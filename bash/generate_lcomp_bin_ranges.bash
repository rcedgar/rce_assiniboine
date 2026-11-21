#!/bin/bash -e

# lcomp = L+fA+fC+fG+fT systematic grid
# old mfe_binning empirical ~equal bins from ROD dataset

mkdir -p ../output
python ../py/generate_lcomp_bin_ranges.py \
	> ../output/lcomp_bin_ranges.tsv
