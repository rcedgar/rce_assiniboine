#!/bin/bash -e

mkdir -p ../big_output
mkdir -p ../tmp_trf

seq 1 64 \
	| parallel -j64 ./run_trf_split.bash

outfa=../big_output/all.trf.fa
rm -f $outfa

cd ../big_splits
rm -f *.dat

for i in `seq 1 64`
do
	echo "masked split$i >> $outfa"
	cat split$i.fa*.mask \
		>> $outfa
done

rm -f *.mask
