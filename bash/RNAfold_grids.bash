#!/bin/bash -e

cat ../output/fasta_grid_names.txt \
	| parallel -j 64 ./RNAfold_grid_split.bash
