#!/bin/bash -e

cd ../big_output

../bin/reorf \
	-getorfs controls.fa \
	-circular \
	-output controls.reorf.fa
