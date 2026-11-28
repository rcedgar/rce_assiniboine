#!/bin/bash -e

cd ../big_rnafold_controls_100_per_family_temps

for T in 0 10 37 60
do
	python /mnt/c/src/py/rnafold_to_label_seq_ss.py T$T.txt \
		> T$T.tsv

	../vleness/vleness \
		--minscore -9 \
		--tsv 1,3 T$T.tsv \
		T$T.vleness.tsv
done

python ../py/joinn.py \
	T0.tsv T10.tsv T37.tsv T60.tsv \
	> joined_T0_T10_T37_T60.tsv


python ../py/joinn.py \
	T0.vleness.tsv T10.vleness.tsv T37.vleness.tsv T37.vleness.tsv \
	> joined_T0_T10_T37_60.vleness.tsv
