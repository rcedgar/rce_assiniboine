#!/bin/bash -e

dir=$1
outputfn=$2

if [ -z "$outputfn" ] ; then
	echo Missing arg >> /dev/stderr
	exit 1
fi

# dir=~/rce_assiniboine/big_clusters/2026-11-06/
# circles.maxlen5000.rotated.pctid90.fa

mkdir -p ../tmp

find $dir | fgrep circles.maxlen5000.rotated.pctid90.fa \
	| sed "-es/\/circles.fa//" \
	> ../tmp/circle_paths.txt

# /watson/rnalab/robedgar/rce_assiniboine/big_clusters/2026-11-06/1112/circles.maxlen5000.rotated.pctid90.fa:8
grep -c "^>" `cat ../tmp/circle_paths.txt` \
	| sed "-es/.*\/\([0-9]*\)\/.*:/\1	/" \
	> $outputfn

ls -lh $outputfn
