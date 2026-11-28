#!/bin/bash -e

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd "$(dirname "$0")/.."
exec /home/bob/miniconda3/envs/null_model/bin/python -u py/avsun_feature_distributions.py
