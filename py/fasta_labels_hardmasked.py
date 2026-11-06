#!/usr/bin/env python3

import sys
import fasta

fafn = sys.argv[1]
max_pct = int(sys.argv[2])

def onseq(label, seq):
	seq = seq.upper()
	L = len(seq)
	if L == 0:
		return
	n = seq.count('N')
	pct = (n*100)/L
	if pct > max_pct:
		print(label)

fasta.ReadSeqsOnSeq(fafn, onseq)
