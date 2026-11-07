#!/usr/bin/python3

import sys
import fev

MAX_EVALUE = 10

FN = sys.argv[1]

#         0     1       2 3
# full_hhr3/24-77/6.8e-07/+
def ParseCmHitStr(Str):
	Fields = Str.split('/')
	if len(Fields) != 4:
		sys.stderr.write("Expected 4 fields got " + Str + "\n")
		assert False

	Name = Fields[0]
	FromTo = Fields[1]
	Evalue = float(Fields[2])
	Strand = Fields[3]
	assert Strand == '+' or Strand == '-'

	Fields2 = FromTo.split('-')
	assert len(Fields2) == 2
	From = int(Fields2[0])
	To = int(Fields2[1])

	Lo = min(From, To)
	Hi = max(From, To)

	Tuple = Name, Lo, Hi, Evalue, Strand
	return Tuple

def GetCmPair():
	cmhits = fev.GetValue("cmhits")
	if cmhits is None:
		return None, None
	Fields = cmhits.split(',')
	assert len(Fields) > 0

	# 3:2:1/full_hhr3/24-77/6.8e-07/+
	FirstField = Fields[0]
	n = FirstField.find('/')
	FirstField = FirstField[:n]
	Fields2 = FirstField.split(':')
	assert len(Fields2) == 3
	N = int(Fields2[0])
	Np = int(Fields2[1])
	Nm = int(Fields2[2])

	assert Fields[0].startswith(FirstField)
	Fields[0] = Fields[0].replace(FirstField, "")[1:]

	CmName_Plus = None
	CmName_Minus = None

	MinEvalue_Plus = 999
	MinEvalue_Minus = 999

	for i in range(0, N):
		Field = Fields[i]
		Name, Lo, Hi, Evalue, Strand = ParseCmHitStr(Field)
		if Evalue > MAX_EVALUE:
			continue
		if Strand == '+':
			if Evalue < MinEvalue_Plus:
				MinEvalue_Plus = Evalue
				CmName_Plus = Name
		elif Strand == '-':
			if Evalue < MinEvalue_Minus:
				MinEvalue_Minus = Evalue
				CmName_Minus = Name
		else:
			assert False

	return CmName_Plus, MinEvalue_Plus, CmName_Minus, MinEvalue_Minus

def OnRec():
	global Len
	global Comment

	Label = fev.Label
	cmname1, evalue1, cmname2, evalue2 = GetCmPair()
	if cmname1 is None or cmname2 is None:
		return
	s = Label
	s += "\tcmplus=%s(%.3g)" % (cmname1, evalue1)
	s += "\tcmminus=%s(%.3g)" % (cmname2, evalue2)
	print(s)

fev.ReadRecsOnRec(FN, OnRec)
