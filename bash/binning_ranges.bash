#!/bin/bash -e

# this is old mfe_binning empirical ~equal bins from ROD dataset
# probably obsolete, see lcomp = L+fA+fC+fG+fT systematic grid

# numpy is in this env; /usr/bin/python does not have it.
py=/home/bob/miniconda3/envs/null_model/bin/python
root="$(cd "$(dirname "$0")/.." && pwd)"

"$py" "$root/py/mfe_binning_ranges.py" \
	--tree "$root/mfe_binning_tree/mfe_binning.tree.tsv" \
	--output "$root/mfe_binning_tree/mfe_binning_ranges.tsv"
