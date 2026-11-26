#!/bin/bash -e

python ../py/rnafold_label_seq_ss_to_circle_features.py \
	controls_folds_100_per_family.tsv \
	> controls_100_per_family_circle_features.tsv
