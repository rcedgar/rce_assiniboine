import sys

# Derived from merge_compos_mfe.py
# Filenames are arguments
# Don't replace "ROD|" -> ""

composfn = sys.argv[1]
mfefn  = sys.argv[2]

label2compos = {}
for line in open(composfn):
	flds = line.rstrip().split('\t')
	label = flds[0]
	label2compos[label] = flds[1:]

found = 0
notfound = 0
f = open(mfefn)
f.readline()
for line in f:
	line = line.rstrip()
	flds = line.split('\t')
	if len(flds) != 2:
		notfound += 1
		continue
	label, MFE = flds
	compos = label2compos.get(label)
	if not compos:
		notfound += 1
		continue
	found += 1
	L, NA, NC, NG, NT = compos
	print(label + "\t" + L + "\t" + MFE + "\t" + NA + "\t" + NC + "\t" + NG + "\t" + NT)

sys.stderr.write("Found %d, not found %d\n" % (found, notfound))
