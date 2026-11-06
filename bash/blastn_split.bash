#!/bin/bash -e

i=$1

if [ -z "$i" ] ; then
	echo "$0: Missing index" >> /dev/stderr
	exit 1
fi

mkdir -p ../big_blast_hits

../bin/blastn \
        -query ../big_splits/split$i.fa \
        -db ../big_output/blastdb/dust_filtered \
        -evalue 1e-6 \
        -outfmt 6 \
        -out ../big_blast_hits/split$i.tsv
