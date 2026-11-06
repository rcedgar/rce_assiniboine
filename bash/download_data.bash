#!/bin/bash -e

mkdir -p ../big_data
cd ../big_data

aws s3 ls s3://assiniboine/fa/control/ \
	| fgrep .fa \
	| sed "-es/.fa//" \
	| sed "-es/.* //" \
	> control_names.txt

aws s3 cp s3://assiniboine/fa/rod/ROD_v1.1.fa.zst .

for name in `cat control_names.txt`
do
	aws s3 cp s3://assiniboine/fa/control/$name.fa .
done

zstd -d ROD_v1.1.fa.zst
rm -f ../big_data/ROD_v1.1.fa.zst

ls -lh ../big_data
