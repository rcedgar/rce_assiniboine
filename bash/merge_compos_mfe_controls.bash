#!/bin/bash -e

cd ../big_output

python ../py/merge_compos_mfe_v2.py \
	controls.compos.tsv \
	controls_label_mfe.tsv \
	> merged_compos_mfe_controls.tsv
