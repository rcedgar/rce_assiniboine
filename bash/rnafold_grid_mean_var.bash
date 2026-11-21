#!/bin/bash -e

mkdir -p ../output
cd ../big_grid_folds

python ../py/rnafold_grid_mean_var.py * \
	> ../output/rnafold_grid_mean_var.tsv
