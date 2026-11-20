#!/bin/bash -e

seq 1 1000 \
	| parallel -j 64 ./RNAfold_control_split.bash
