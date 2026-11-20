#!/bin/bash -e

infa=dust_masked_gt25pct.fa

cd ../big_output

../bin/reorf \
	-getorfs $infa \
	-circular \
	-longest_orf \
	-sep : \
	-output dust_masked_gt25pct.reorf_longest.fa

usearch \
	-sortbylength dust_masked_gt25pct.reorf_longest.fa \
	-minseqlength 30 \
	-fastaout dust_masked_gt25pct.reorf_longest.bylength.fa

usearch \
	-fastx_uniques dust_masked_gt25pct.reorf_longest.bylength.fa \
	-fastaout dust_masked_gt25pct.reorf_longest.bylength.uniques.fa \
	-uc dust_masked_gt25pct.reorf_longest.bylength.uniques.uc

usearch \
	-sortbylength dust_masked_gt25pct.reorf_longest.bylength.uniques.fa \
	-fastaout dust_masked_gt25pct.reorf_longest.bylength.uniques.bylength.fa
