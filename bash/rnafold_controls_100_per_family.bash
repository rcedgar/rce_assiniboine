#!/bin/bash -e

outdir=../big_rnafold_controls_100_per_family_temps/
mkdir -p $outdir
cd $outdir

RNAfold \
	--infile=../big_output/controls_100_per_family.fa \
	--circ \
	--noPS \
	-T 0 \
	--jobs=32 \
	--outfile=T0.txt \
	> T0.stdout \
	2> T0.stderr

RNAfold \
	--infile=../big_output/controls_100_per_family.fa \
	--circ \
	--noPS \
	-T 10 \
	--jobs=32 \
	--outfile=T10.txt \
	> T10.stdout \
	2> T10.stderr

RNAfold \
	--infile=../big_output/controls_100_per_family.fa \
	--circ \
	--noPS \
	--jobs=32 \
	--outfile=T37.txt \
	> T37.stdout \
	2> T37.stderr

RNAfold \
	--infile=../big_output/controls_100_per_family.fa \
	--circ \
	--noPS \
	--jobs=32 \
	--outfile=p.txt \
	-p \
	> p.stdout \
	2> p.stderr
