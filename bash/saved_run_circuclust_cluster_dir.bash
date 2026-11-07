#!/bin/bash -e

dir=~/rce_assiniboine/big_clusters/2026-11-06/

cd $dir
find | fgrep circles.fa \
	| sed "-es/\/circles.fa//" \
	> circle_paths.txt

for path in `cat circle_paths.txt`
do
	echo
	echo $path
	echo
	cd $path
	usearch \
		-sortbylength circles.fa \
		-maxseqlength 5000 \
		-fastaout circles.maxlen5000.fa

	n=`cat circles.maxlen5000.fa | grep "^>" | wc -l`
	if [ $n -gt 1 ] ; then
		circuclust \
			-cluster circles.maxlen5000.fa \
			-id 0.9 \
			-fastaout circles.maxlen5000.rotated.pctid90.fa \
			-tsvout circuclust.tsv
	fi
	cd $dir
done
