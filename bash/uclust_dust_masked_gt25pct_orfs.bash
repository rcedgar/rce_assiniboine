#!/bin/bash -e

cd ../big_output

infa=dust_masked_gt25pct.reorf_longest.bylength.uniques.bylength.fa

usearch \
	-cluster_fast $infa \
	-id 0.7 \
	-sizeout \
	-centroids orfs_dust_masked_gt25pct.id70.fa \
	-uc orfs_dust_masked_gt25pct.id70.uc
