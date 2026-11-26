#!/bin/bash -e

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd ../big_output

mkdir -p ../analyze_circle_features_filtered

python ../py/analyze_circle_features.py \
	controls_circle_features_filtered.tsv \
	../analyze_circle_features_filtered
