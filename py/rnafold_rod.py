#!/usr/bin/python3

import sys

# Input is tsv with i,j pairs only; first line is sequence.
# To make input:
#	rnafold_pairprobs.py sequence_0001_dp.ps > pairprobs.tsv

FileName = sys.argv[1]

Pairs = []
f = open(FileName)
Seq = f.readline().strip()
for c in Seq:
	if not c in "ACGTUN":
		sys.stderr.write("Invalid letter '%c' in RNA sequence\n" % c)
L = len(Seq)

while 1:
	Line = f.readline()
	if len(Line) == 0:
		break
	Fields = Line.strip().split()
	assert len(Fields) >= 2

	i = int(Fields[0])
	j = int(Fields[1])
	P = float(Fields[2])
	if P < 0.5:
		continue
	assert i >= 0 and i < L
	assert j >= 0 and j < L
	Pairs.append((i, j))

ItoJ = [ None ] * (L+1)

for Pair in Pairs:
	i, j = Pair
	assert i <= L and j <= L
	if ItoJ[i] == None:
		ItoJ[i] = j
	else:
		assert ItoJ[i] == j

'''
         0                               L-1                           2*L-1
		 |                                |                              |
         ==================================+++++++++++++++++++++++++++++++
		         |               |         |         |
		         i               j         L        L+i
				 <--- (j - i) --->
				                 <--- (L + i - j) --->
'''
def GetDist12(i, j):
	assert i < j
	d1 = j - i
	assert d1 >= 0 and d1 < L

	d2 = L + i - j
	assert d2 >= 0 and d2 < L

	return min(d1, d2)

def GetDist(i, j):
	if i < j:
		return GetDist12(i, j)
	elif i > j:
		return GetDist12(j, i)
	else:
		return 0

Maxd = 0
Maxd_i = 0
Maxd_j = 0
Ds = L*[ 0 ]
for i in range(0, L):
	j = ItoJ[i]
	if j == None:
		continue
	d = GetDist(i, j)
	Ds[i] = d
	if d > Maxd:
		Maxd = d
		Max_i = i
		Max_j = j

if 1:
	for i in range(0, L):
		d = Ds[i]
		print(("d\t%d\t%d" % (i, d)))

Equator1 = min(Max_i, Max_j)
Equator2 = max(Max_i, Max_j)
Pole1 = (Equator1 + Equator2)/2
Pole2 = (Pole1 + L//2)%L

M = 0
Good = 0
for i in range(0, L):
	j = ItoJ[i]
	if j == None:
		continue
	M += 1

	Actual_d = GetDist(i, j)

	DistPole1_i = GetDist(i, Pole1)
	DistPole2_i = GetDist(i, Pole2)
	MinDistPole_i = min(DistPole1_i, DistPole2_i)

	Perfect_d = MinDistPole_i*2
	Error = abs(Actual_d - Perfect_d)
	if Error < L/10:
		Good += 1

#	print("i=%d j=%d d=%d perfd=%d" % (i, j, Actual_d, Perfect_d))

Score = Good/M

print(("Rodness=%.4f	L=%d Equator=%d,%d	Poles=%d,%d	File=%s" % (Score, L, Equator1, Equator2, Pole1, Pole2, FileName)))
