#!/bin/bash -e

# Input is .fold file from RNAfold, format is fasta-like with
# sequence label followed by sequence followed by one or more lines
# of ss in ((.)) format.
# dot-bracket line ends with MFE in parentheses e.g. ))))))))) (-25.50)

# Output is TSV with 1=label 2=MFE

import sys

label = None
for line in open(sys.argv[1]):
	if line.startswith('>'):
		label = line[1:-1]
		continue
	if line[-2] == ")":
		# print(line)
		line = line.strip()
		flds = line.split()
		MFE = flds[1][1:-1]
		print(label + "\t" + MFE)
		label = None
		continue
