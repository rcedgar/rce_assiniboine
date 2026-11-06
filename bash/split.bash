#!/bin/bash -e

N=64

mkdir -p ../big_splits
cd ../big_splits

../bin/usearch \
	-fastx_split ../big_output/dust_filtered.fa \
	-splits $N \
	-outname split@.fa
