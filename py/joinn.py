#!/usr/bin/python3

# join lines from two or more tsv files
# primary key is first field
# no assumption about sort
# assumes key found once per file, dupes overwrite
# include keys found only in one file
# missing field is .

import sys

NA = "."

FNs = sys.argv[1:]

FieldCounts = []

FirstKey = None

Keys = set()
def ReadFile(FN):
	global FirstKey
	KeyToFields = {}
	FieldCount = None
	for Line in open(FN):
		Fields = Line[:-1].split('\t')
		Key = Fields[0]
		if FirstKey is None:
			FirstKey = Key
		Keys.add(Key)
		Fields = Fields[1:]
		if FieldCount is None:
			FieldCount = len(Fields)
		else:
			assert len(Fields) == FieldCount
		KeyToFields[Key] = Fields
	return FieldCount, KeyToFields

KeyToFieldsVec = []
for FN in FNs:
	FieldCount, KeyToFields = ReadFile(FN)
	FieldCounts.append(FieldCount)
	KeyToFieldsVec.append(KeyToFields)

KeyList = list(Keys)
KeyList.sort()

N = len(KeyToFieldsVec)
assert len(FieldCounts) == N
for i in range(len(KeyList)):
	if i == 0:
		Key = FirstKey
	else:
		Key = KeyList[i]
		if Key == FirstKey:
			continue
	s = Key
	for i in range(N):
		FieldCount = FieldCounts[i]
		KeyToFields = KeyToFieldsVec[i]
		try:
			Fields = KeyToFields[Key]
		except:
			Fields = None
		if Fields is None:
			Fields = []
			for j in range(FieldCount):
				Fields.append(NA)

		for Field in Fields:
			s += "\t" + Field
	print(s)
