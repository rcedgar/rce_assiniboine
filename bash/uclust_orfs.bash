#!/bin/bash -e

cd ../big_output

usearch \
	-cluster_fast dust_filtered.reorf.fa \
	-id 0.7 \
	-sizeout \
	-centroids orfs.id70.fa \
	-uc orfs.id70.uc
