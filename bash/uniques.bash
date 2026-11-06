#!/bin/bash -e

cd ../big_output

usearch \
	-fastx_uniques all.fa \
	-fastaout uniques.fa \
	-uc uniques.uc
