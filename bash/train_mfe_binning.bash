#!/bin/bash -e

# numpy is in this env; /usr/bin/python does not have it.
py=/home/bob/miniconda3/envs/null_model/bin/python
root="$(cd "$(dirname "$0")/.." && pwd)"

cd "$root/big_output"

"$py" "$root/py/mfe_binning_build_bins.py" \
	--input merged_compos_mfe.tsv \
	--out-prefix mfe_binning \
	--n-min 2000 \
	--tau 0.02 \
	--max-depth 40 \
	--holdout-frac 0.2 \
	--seed 1
