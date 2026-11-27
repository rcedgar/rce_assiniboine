#!/bin/bash -e

cd "$(dirname "$0")/.."

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PY=/home/bob/miniconda3/envs/null_model/bin/python

"$PY" -c "import sys; sys.path.insert(0, 'py'); import dotbracket_circular_features as m; m._check_rod_features()"

"$PY" -u py/rnafold_label_seq_ss_to_circle_features.py \
	big_output/controls_label_seq_ss.tsv \
	> big_output/controls_circle_features.tsv

(
	cd py
	"$PY" -u filter_circle_features_VLE_outliers.py
)

"$PY" -u py/analyze_circle_features.py \
	big_output/controls_circle_features_filtered.tsv \
	analyze_circle_features_filtered

"$PY" -u py/construct_VLE_discriminator_score.py \
	big_output/controls_circle_features_filtered.tsv \
	construct_VLE_discriminator_score
