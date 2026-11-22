#!/bin/bash -e

cd "$(dirname "$0")/.."

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec /home/bob/miniconda3/envs/null_model/bin/python -u py/compare_z_mfe_bins.py
