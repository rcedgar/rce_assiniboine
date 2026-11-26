#!/bin/bash -e

cd ../big_output

python ../py/rnafold_label_seq_ss_to_circle_features.py \
	controls_label_seq_ss.tsv \
	> controls_circle_features.tsv
