#!/bin/bash -e

cd "$(dirname "$0")/.."
python3 py/select_control_labels.py
