#!/bin/bash -e

cd "$(dirname "$0")/.."

exec /home/bob/miniconda3/envs/null_model/bin/python -u py/scatter_rod_subset_z.py
