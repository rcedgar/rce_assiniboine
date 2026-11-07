#!/bin/bash -e

dir=~/rce_assiniboine/big_clusters/2026-11-06/

mkdir -p ../tmp

find $dir | fgrep circles.fa \
	| sed "-es/\/circles.fa//" \
	> ../tmp/circle_paths.txt

cat ../tmp/circle_paths.txt \
	| parallel -j 10 \
	./run_cmscan_cluster_path.bash
