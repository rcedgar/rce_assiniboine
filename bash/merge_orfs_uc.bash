#!/bin/bash -e

cd ../big_output

python ../py/uc_merge2.py dust_filtered.reorf_longest.bylength.uniques.uc orfs.id70.uc \
	> orfs.id70.merged.tsv
