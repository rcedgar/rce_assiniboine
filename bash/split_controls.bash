#!/bin/bash -e

rm -rf ../big_control_fas_split
mkdir -p ../big_control_fas_split

cd ../big_output

python $py/fasta_shuffle_order.py controls_length100to6000_bylength.fa \
	> controls_length100to6000_shuffled.fa

usearch \
    -fastx_split ../big_output/controls_length100to6000_shuffled.fa \
    -splits 1000 \
    -outname ../big_control_fas_split/split@.fa
