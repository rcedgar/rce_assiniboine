diamond \
	blastx \
	--sensitive \
	--threads 64 \
	-d centroid37781 \
	-q ../big_data/LSC_v1.1.fa \
	-o hits.tsv
