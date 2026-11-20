#!/bin/bash -e

n=$1

fa=../big_control_fas_split/split$n.fa

if [ ! -s $fa ] ; then
	echo "Not found fa=$fa"
	exit 1
fi

RNAfold --circ $fa > ../big_control_folds/split$n.txt
