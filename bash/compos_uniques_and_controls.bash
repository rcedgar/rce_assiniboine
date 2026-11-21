#!/bin/bash -e

cd ../big_output

reorf -compos uniques.fa -output uniques.compos.tsv

reorf -compos controls.fa -output controls.compos.tsv
