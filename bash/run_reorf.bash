#!/bin/bash -e

cd ../big_output

../bin/reorf \
	-getorfs dust_filtered.fa \
	-circular \
	-output dust_filtered.reorf.fa
