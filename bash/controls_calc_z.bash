#!/bin/bash -e

cd ../big_output

python ../py/mfe_binning_calc_z.py \
	--tree ../mfe_binning_tree/mfe_binning.tree.tsv \
	--input merged_compos_mfe_controls.tsv \
	| sort -gk5 \
	> controls_calc_z.tsv
