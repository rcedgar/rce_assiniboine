#!/bin/bash -e

echo ##################
echo "OBSOLETE I THINK"
echo ##################
exit 1

python ../py/rnafold_to_label_seq_ss.py \
	control_RNAfold.txt  \
	> control_label_seq_ss.tsv
