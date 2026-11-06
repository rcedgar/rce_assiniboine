#!/bin/bash -e

cd ../big_output

infa=dust_filtered.reorf_longest.bylength.uniques.bylength.fa

usearch \
	-cluster_fast $infa \
	-id 0.7 \
	-sizeout \
	-centroids orfs.id70.fa \
	-uc orfs.id70.uc
