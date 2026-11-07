#!/bin/bash -e

anadir=~/rce_assiniboine/analysis/2026-11-06

./species_counts_cluster_dir.bash \
	~/rce_assiniboine/big_clusters/2026-11-06/ \
	~/rce_assiniboine/analysis/2026-11-06/cluster_speciescount.tsv

python ../py/fev_append_species_counts.py \
	$anadir/top200_selected.fev \
	$anadir/cluster_speciescount.tsv \
	> $anadir/top200_selected_with_species_count.fev

ls -lh $anadir/top200_selected_with_species_count.fev
