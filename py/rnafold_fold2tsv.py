#!/usr/bin/python3

import sys
import math
import fasta

# Input is .fold file from RNAfold, format is fasta-like with
# sequence label followed by sequence followed by one or more lines
# of ss in ((.)) format.

# Output is tsv with i,j pairs

FileName = sys.argv[1]

def SetPairs():
	global L, Pairs, ItoJ, Parens

	Pairs = []
	L = len(Parens)
	Stack = []
	for i in range(0, L):
		c = Parens[i]
		if c == '(':
			Stack.append(i)
		elif c == ')':
			i2 = Stack.pop()
			Pair = (i2, i)
			Pairs.append(Pair)
	assert len(Stack) == 0

def GetDist12(i, j):
	global L, Pairs, ItoJ, Parens

	assert i < j
	d1 = j - i
	assert d1 >= 0 and d1 < L

	d2 = L + i - j
	assert d2 >= 0 and d2 < L

	return min(d1, d2)

def GetDist(i, j):
	global L, Pairs, ItoJ, Parens

	if i < j:
		return GetDist12(i, j)
	elif i > j:
		return GetDist12(j, i)
	else:
		return 0

def SetItoJ():
	global L, Pairs, ItoJ

	ItoJ = [ None ] * (L+1)

	for Pair in Pairs:
		i, j = Pair
		assert i <= L and j <= L
		if ItoJ[i] == None:
			ItoJ[i] = j
		else:
			assert ItoJ[i] == j

TWOPI = 2.0*3.1415927

def Output():
	global L, Pairs, ItoJ, Parens, Label
	global Label, Seq

	for i in range(0, L):
		Theta = TWOPI*float(i)/float(L)
		x = (1.0 + math.cos(Theta))/2.0
		y = (1.0 + math.sin(Theta))/2.0
		j = ItoJ[i]
		if j == None:
			S = 1
			s = "%d" % i
			s += "\t%.4f" % x
			s += "\t%.4f" % y
			s += "\t."
			s += "\t."
			print(s)
			continue

		S = 1
		s = "%d" % i
		s += "\t%.4f" % x
		s += "\t%.4f" % y
		s += "\t%.4f" % S
		s += "\t%d" % j
		print(s)

FirstLine = False
SecondLine = False
print("i\tx\ty\tS\tj")
for Line in open(FileName):
	if Line.startswith(">"):
		Label = Line[1:-1]
		if Label.endswith("_0001"):
			Label = Label.replace("_0001", "")
		FirstLine = True
		continue
	elif FirstLine:
		Seq = Line[:-1]
		SecondLine = True
		FirstLine = False
	elif SecondLine:
		Parens = Line[:-1].split()[0]
		L = len(Parens)
		SetPairs()
		SetItoJ()
		Output()
		sys.exit(0)
assert False
