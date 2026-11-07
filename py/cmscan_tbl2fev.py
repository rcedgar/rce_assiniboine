#!/usr/bin/python3

import sys

'''
#target name         accession query name                            accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#------------------- --------- ------------------------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
#                  0         1                                     2         3   4        5        6        7        8      9    10   11   12    13     14        15  16                    17
hhrbz_dv4            -         NODE_67741_length_1127_cov_9.528463 -          cm        1       70      632      579      -    no    1 0.48   0.0   26.7   7.5e-07 !   -
hhrbz_dv4            -         NODE_67741_length_1127_cov_9.528463 -          cm        1       70      214      286      +    no    1 0.51   0.0   22.0   1.9e-05 !   -
'''

TblFN = sys.argv[1]
MIN_EVALUE = float(sys.argv[2])

def argsort(arr, reverse=False):
    return sorted(range(len(arr)), key=lambda i: arr[i], reverse=reverse)

HitTuples = []
CurrentLabel = None

def SortHitTuples(HitTuples):
	Order = []
	SelectedHitTuples = []
	N = len(HitTuples)
	if N == 0:
		return 0, 0, Order, SelectedHitTuples
	
	Np = 0
	Nm = 0
	IndexToEvalue = []
	Index = 0
	for i in range(0, N):
		Tuple = HitTuples[i]
		Model, SeqLo, SeqHi, HitStrand, Evalue = Tuple
		if Evalue <= MIN_EVALUE:
			IndexToEvalue.append(Evalue)
			SelectedHitTuples.append(Tuple)
			Index += 1
			if HitStrand == '+':
				Np += 1
			elif HitStrand == '-':
				Nm += 1
			else:
				assert False
	Order = argsort(IndexToEvalue, reverse=True)
	return Np, Nm, Order, SelectedHitTuples

def GetHitStr(HitIndex):
	global HitTuples
	if HitIndex == None:
		return None
	Model, SeqLo, SeqHi, HitStrand, Evalue = HitTuples[HitIndex]

	HitStr = "%s/%d-%d/%.3g/%s" % (Model, int(SeqLo), int(SeqHi), Evalue, HitStrand)
	return HitStr

def Flush():
	global CurrentLabel, HitTuples
	
	N = len(HitTuples)
	if N == 0:
		return

	Np, Nm, Order, SelectedHitTuples = SortHitTuples(HitTuples)
	M = len(Order)
	assert len(SelectedHitTuples) == M
	assert Np + Nm == M

	Str = ""
	for k in range(0, M):
		HitIndex = Order[k]
		if k > 0:
			Str += ","
		Str += GetHitStr(HitIndex)

	s = "label=" + CurrentLabel
	s += "\tcmhits=%d:%d:%d/" % (N, Np, Nm)
	s += Str
	print(s)

for Line in open(TblFN):
	if Line.startswith('#'):
		continue
	Fields = Line.split()
	if len(Fields) < 18:
		sys.stderr.write("WARNING < 18 fields in %s\n" % TblFN)
		continue
	QueryLabel = Fields[2]
	Model = Fields[0]
	SeqLo = Fields[7]
	SeqHi = Fields[8]
	Strand = Fields[9]
	Evalue = float(Fields[15])
	if Evalue > MIN_EVALUE:
		continue

	HitTuple = (Model, SeqLo, SeqHi, Strand, Evalue)

	if QueryLabel != CurrentLabel:
		Flush()
		CurrentLabel = QueryLabel
		HitTuples = []

	HitTuples.append(HitTuple)

Flush()
