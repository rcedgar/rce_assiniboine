#!/bin/bash -e

cd "$(dirname "$0")/.."

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec /home/bob/miniconda3/envs/null_model/bin/python -u py/predict_z_7cols.py \
	--input big_output/merged_compos_mfe_7cols.head1000.tsv \
	--output big_output/merged_compos_mfe_7cols.head1000_predicted_z.tsv
