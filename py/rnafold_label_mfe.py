#!/bin/bash -e

# Input is .fold file from RNAfold, format is fasta-like with
# sequence label followed by sequence followed by one or more lines
# of ss in ((.)) format.
# dot-bracket line ends with MFE in parentheses e.g. ))))))))) (-25.50)

# Output is TSV with 1=label 2=MFE

import sys
import re

pattern = r"\(\s*([+-]?(?:\d+\.?\d*|\.\d+))\s*\)"

for fn in sys.argv[1:]:
	label = None
	for line in open(fn):
		if line.startswith('>'):
			label = line[1:-1]
			continue
		if line[-2] == ")":
			line = line.strip()
			match = re.search(pattern, line)
			if not match:
				continue
			MFE = match.group(1)
			print(label + "\t" + MFE)
			label = None
			continue
