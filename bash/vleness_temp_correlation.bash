#!/bin/bash -e

cd "$(dirname "$0")/.."

exec /home/bob/miniconda3/envs/null_model/bin/python -u py/vleness_temp_correlation.py
