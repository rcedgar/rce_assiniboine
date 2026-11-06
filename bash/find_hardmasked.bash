#!/bin/bash -e

cd ../big_output

pct=25

python ../py/fasta_labels_hardmasked.py uniques.dust.fa $pct \
	> dust_masked_gt${pct}pct.labels
