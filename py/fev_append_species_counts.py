#!/usr/bin/env python3

import sys

fevfn=sys.argv[1]
countsfn=sys.argv[2]

cluster2n = {}
for line in open(countsfn):
	flds = line[:-1].split('\t')
	assert len(flds) == 2
	cluster = flds[0]
	n = flds[1]
	cluster2n[cluster] = n

for line in open(fevfn):
	line = line[:-1]
	flds = line.split('\t')
	assert flds[0].startswith("cluster_nr=")
	cluster = flds[0].split('=')[1]
	n = cluster2n.get(cluster, "0")
	if n != "0":
		line += "\tspecies=" + n
	print(line)
