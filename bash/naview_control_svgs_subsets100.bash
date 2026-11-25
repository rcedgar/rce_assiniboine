#!/bin/bash -e

cd ../big_output

outdir=../big_control_svgs_100_per_family
outdir2=../big_control_svgs_100_per_family_naview_only

mkdir -p $outdir $outdir2

../build-linux/rce_naview \
	--circular \
	--naview-only \
	--rainbow \
	--threads 32 \
	--out-dir $outdir2 \
	controls_folds_100_per_family.tsv

../build-linux/rce_naview \
	--circular \
	--rainbow \
	--threads 32 \
	--out-dir $outdir \
	controls_folds_100_per_family.tsv
