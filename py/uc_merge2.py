#!/usr/bin/env python3

# Merge two .uc files from usearch cluster*/uniques
# output is TSV 1=centroid 2=member 3=pctid

import sys

ucfn1 = sys.argv[1]
ucfn2 = sys.argv[2]

centroid2members1 = {}
for line in open(ucfn1):
	flds = line[:-1].split('\t')
	X = flds[0] # S/H/C...
	if X == 'S':
		q = flds[8]
		if q == "*":
			print(line)
			assert False
		centroid2members1[q] = []
		continue
	elif X == 'H':
		q = flds[8]
		t = flds[9]
		if q == "*" or t == "*":
			print(line)
			assert False
		centroid2members1[t].append(q)

def output_members(t, pctid):
	members = centroid2members1.get(t)
	if members:
		for q in members:
			if q == "*" or t == "*":
				print(line)
				assert False
			print("H\t" + q + "\t" + t + "\t" + pctid)

for line in open(ucfn2):
	flds = line[:-1].split('\t')
	X = flds[0]
	if X == 'S':
		q = flds[8]
		if q == "*":
			print(line)
			assert False
		print("S\t" + q + "\t" + q + "\t100")
		output_members(q, "100.00")
	elif X == 'H':
		pctid = flds[3]
		q = flds[8]
		t = flds[9]
		if q == "*" or t == "*":
			print(line)
			assert False
		print("H\t" + q + "\t" + t + "\t" + pctid)
		output_members(q, pctid)
