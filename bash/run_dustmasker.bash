#!/bin/bash -e

cd ../big_output

../bin/dustmasker \
	-in uniques.fa \
	-hard_masking \
	-outfmt fasta \
	-out uniques.dust.fa

ls -lh uniques.dust.fa
