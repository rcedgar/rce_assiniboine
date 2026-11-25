#!/bin/bash -e

cd "$(dirname "$0")/.."
python3 py/extract_control_folds.py
