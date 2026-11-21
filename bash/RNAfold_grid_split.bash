#!/bin/bash -e

name=$1

fa=../big_fasta_grid/$name.fa

if [ ! -s $fa ] ; then
	echo "Not found fa=$fa"
	exit 1
fi

RNAfold --circ --noPS $fa > ../big_grid_folds/$name.txt
