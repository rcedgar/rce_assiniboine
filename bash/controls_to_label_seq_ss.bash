#!/bin/bash -e

cd ../big_output

python ../py/rnafold_to_label_seq_ss.py \
	control_RNAfold.txt \
	> controls_label_seq_ss.tsv
