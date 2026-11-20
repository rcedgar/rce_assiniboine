#!/usr/bin/python3

import sys

# Input is tsv with optional header:
#	sample_name	length	mfe	dot_bracket
# Dot-bracket is converted to 0-based i,j pairs as in rnafold_fold2tsv.py.
# Rodness is the score from rnafold_rod.py.

def GetDist12(i, j):
	global L
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

def PairsFromDotBracket(Parens):
	Pairs = []
	Stack = []
	n = len(Parens)
	for i in range(0, n):
		c = Parens[i]
		if c == '(':
			Stack.append(i)
		elif c == ')':
			assert len(Stack) > 0
			i2 = Stack.pop()
			Pairs.append((i2, i))
		elif c == '.':
			continue
		else:
			assert False
	assert len(Stack) == 0
	return Pairs

def Rodness(Parens):
	global L

	L = len(Parens)
	Pairs = PairsFromDotBracket(Parens)

	ItoJ = [ None ] * (L+1)
	for Pair in Pairs:
		i, j = Pair
		assert i <= L and j <= L
		if ItoJ[i] == None:
			ItoJ[i] = j
		else:
			assert ItoJ[i] == j

	Maxd = 0
	Max_i = None
	Max_j = None
	for i in range(0, L):
		j = ItoJ[i]
		if j == None:
			continue
		d = GetDist(i, j)
		if d > Maxd:
			Maxd = d
			Max_i = i
			Max_j = j

	assert Max_i != None and Max_j != None

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

	assert M > 0
	Score = Good/M
	return Score, Equator1, Equator2, Pole1, Pole2

def Main():
	FileName = sys.argv[1]

	print("sample_name\tlength\tmfe\trodness\tequator1\tequator2\tpole1\tpole2")

	SeenData = False
	for Line in open(FileName):
		Line = Line.strip()
		if len(Line) == 0:
			continue
		Fields = Line.split('\t')
		if not SeenData:
			try:
				int(Fields[1])
			except ValueError:
				continue
		SeenData = True

		assert len(Fields) == 4
		Name = Fields[0]
		Length = int(Fields[1])
		Mfe = Fields[2]
		Parens = Fields[3]
		assert Length == len(Parens)

		Score, Equator1, Equator2, Pole1, Pole2 = Rodness(Parens)
		print("%s\t%d\t%s\t%.4f\t%d\t%d\t%d\t%d" %
			(Name, Length, Mfe, Score, Equator1, Equator2, Pole1, Pole2))

if __name__ == "__main__":
	Main()
