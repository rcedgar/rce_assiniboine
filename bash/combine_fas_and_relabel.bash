#!/bin/bash -e

cd ../big_data

mkdir -p ../big_output

cat ROD_v1.1.fa \
	| sed "-es/ .*//" \
	| sed "-es/^>/>ROD|/" \
	> ../big_output/all.fa

for name in `cat control_names.txt`
do
	echo "=== $name ==="
	cat $name.fa \
		| sed "-es/ .*//" \
		| sed "-es/^>/>$name|/" \
		>> ../big_output/all.fa
done
