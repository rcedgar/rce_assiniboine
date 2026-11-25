#!/bin/bash -e

# Input is .fold files from RNAfold, format is fasta-like with
# sequence label followed by sequence followed by one or more lines
# of ss in ((.)) format.
# dot-bracket line ends with MFE in parentheses e.g. ))))))))) (-25.50)

# Output is TSV with 1=label, 2=seq (T->U), 3=dot-bracket ss

import sys
import re

for fn in sys.argv[1:]:
	label = None
	for line in open(fn):
		if line.startswith('>'):
			if not label is None:
				print(label + "\t" + seq + "\t" + ss)
			label = line[1:-1]
			continue
		if line[-2] == ")":
			ss = line.strip().split()[0] # remove MFE
			continue
		if line[0].isupper():
			seq = line.strip()
if not label is None:
	print(label + "\t" + seq + "\t" + ss)
