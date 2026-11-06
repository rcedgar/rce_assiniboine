#!/bin/bash -e

cd ../big_output

../bin/blastn \
        -query dust_filtered.fa \
        -db blastdb/dust_filtered \
        -evalue 1e-6 \
        -num_threads 32 \
        -outfmt 6 \
        -out blast_hits.tsv
