#!/bin/bash -e

dir=../big_clusters/2026-11-06/
rm -rf $dir

python ../py/make_cluster_dir.py \
	../analysis/2026-11-06/top200_selected.fev \
	../big_output/orfs.id70.merged.tsv \
	../big_output/all.fa \
	../big_output/dust_filtered.reorf_longest.fa \
	$dir
