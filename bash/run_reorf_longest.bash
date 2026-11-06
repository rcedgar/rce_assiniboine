#!/bin/bash -e

cd ../big_output

../bin/reorf \
	-getorfs dust_filtered.fa \
	-circular \
	-longest_orf \
	-sep : \
	-output dust_filtered.reorf_longest.fa

usearch \
	-sortbylength dust_filtered.reorf_longest.fa \
	-minseqlength 30 \
	-fastaout dust_filtered.reorf_longest.bylength.fa

usearch \
	-fastx_uniques dust_filtered.reorf_longest.bylength.fa \
	-fastaout dust_filtered.reorf_longest.bylength.uniques.fa \
	-uc dust_filtered.reorf_longest.bylength.uniques.uc

usearch \
	-sortbylength dust_filtered.reorf_longest.bylength.uniques.fa \
	-fastaout dust_filtered.reorf_longest.bylength.uniques.bylength.fa
