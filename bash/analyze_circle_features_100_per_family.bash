#!/bin/bash -e

cd "$(dirname "$0")/.."

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec /home/bob/miniconda3/envs/null_model/bin/python -u py/analyze_circle_features.py \
	big_output/controls_folds_100_per_family_circle_features_VLE_excludes_discarded.tsv \
	analyze_circle_features_100_per_family
