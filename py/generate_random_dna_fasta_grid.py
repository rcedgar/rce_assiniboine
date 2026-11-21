#!/usr/bin/env python3

import sys
import os

'''
# head ../output/lcomp_bin_ranges.tsv
# ../py/generate_lcomp_bin_ranges.py
# L     fA      fC      fG      fT
240     0.000   0.000   0.125   0.875
240     0.000   0.000   0.250   0.750
240     0.000   0.000   0.375   0.625
240     0.000   0.000   0.500   0.500
240     0.000   0.000   0.625   0.375
'''

outdir = "../big_fasta_grid/"
rangefn = "../output/lcomp_bin_ranges.tsv"
seqs_per_range = 100
pyfn = "../py/generate_random_dna_fasta.py"
	# parser.add_argument("label", help="FASTA defline label")
	# parser.add_argument("L", type=int, help="Sequence length")
	# parser.add_argument("fA", type=float, help="Frequency of A")
	# parser.add_argument("fC", type=float, help="Frequency of C")
	# parser.add_argument("fG", type=float, help="Frequency of G")
	# parser.add_argument("fT", type=float, help="Frequency of T")
	# parser.add_argument("seed", type=float, help="Random number seed (integer)")


rc = os.system("mkdir -p " + outdir)
assert rc==0

for line in open(rangefn):
	if line.startswith('#'):
		continue
	L, fA, fC, fG, fT = line.rstrip().split('\t')
	rangename = "L" + L + "_" + fA + "_" + fC + "_" + fG + "_" + fT
	fafn = outdir + rangename
	print(fafn)
	rc = os.system("rm -f " + fafn)
	for seq_idx in range(seqs_per_range):
		label = rangename + "_seed" + str(seq_idx)
		cmdline = pyfn
		cmdline += " " + label
		cmdline += " " + L
		cmdline += " " + fA
		cmdline += " " + fC
		cmdline += " " + fG
		cmdline += " " + fT
		cmdline += " " + str(seq_idx) # seed
		cmdline += " >> " + fafn
		rc = os.system(cmdline)
		assert rc == 0