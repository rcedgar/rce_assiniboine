#!/bin/bash -e

cd ../big_output

python ../py/rnafold_label_mfe.py control_RNAfold.txt \
	> controls_label_mfe.tsv
