#!/bin/bash -e

cd ../big_output

mkdir -p blastdb

makeblastdb \
	-in dust_filtered.fa \
	-dbtype nucl \
	-out blastdb/dust_filtered

ls -lh blastdb
