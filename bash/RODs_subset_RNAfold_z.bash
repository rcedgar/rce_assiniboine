#!/bin/bash -e

cd "$(dirname "$0")/.."

export PATH="/home/bob/miniconda3/envs/std/bin:$PATH"

exec /home/bob/miniconda3/envs/null_model/bin/python -u py/rnafold_fasta_z.py \
	--input big_output/RODs_subset.fa \
	--output big_output/RODs_subset_RNAfold_z.tsv
