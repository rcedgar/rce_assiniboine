#!/bin/bash -e

# numpy and pytest are in this env; /usr/bin/python does not have them.
py=/home/bob/miniconda3/envs/null_model/bin/python

cd "$(dirname "$0")/.."
"$py" -m pytest py/mfe_binning_test.py
