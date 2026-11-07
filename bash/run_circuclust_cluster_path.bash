#!/bin/bash -e

path=$1

echo $path

usearch \
	-sortbylength $path/circles.fa \
	-maxseqlength 5000 \
	-quiet \
	-fastaout $path/circles.maxlen5000.fa

n=`cat $path/circles.maxlen5000.fa | grep "^>" | wc -l`
if [ $n -gt 1 ] ; then
	circuclust \
		-cluster $path/circles.maxlen5000.fa \
		-quiet \
		-id 0.9 \
		-fastaout $path/circles.maxlen5000.rotated.pctid90.fa \
		-tsvout $path/circuclust.tsv
fi
