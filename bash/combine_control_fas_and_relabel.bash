#!/bin/bash -e

mkdir -p ../output
mkdir -p ../big_output

outfa=../big_output/controls.fa

rm -f $outfa

for name in `cat ../rce_edit/VLE_control_names.txt`
do
	echo "=== VLE $name ==="
	cat ../big_data/$name.fa \
		| sed "-es/^>/>VLE_${name}:/" \
		| sed "-es/>\([^ ]*\) .*/>\1/" \
		>> $outfa
done

for name in `cat ../rce_edit/DECOY_control_names.txt`
do
	echo "=== DECOY $name ==="
	cat ../big_data/$name.fa \
		| sed "-es/^>/>DECOY_${name}:/" \
		| sed "-es/>\([^ ]*\) .*/>\1/" \
		>> $outfa
done
