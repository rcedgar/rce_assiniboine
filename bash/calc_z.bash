#!/bin/bash -e

python ../py/mfe_binning_calc_z.py \
	--tree ../mfe_binning_tree/mfe_binning.tree.tsv \
	--input merged_compos_mfe_7cols.tsv \
	| sort -gk5 \
	> calc_z.tsv
