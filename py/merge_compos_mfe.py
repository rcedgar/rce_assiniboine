import sys

'''
# cut -f1,2,3 ~/assiniboine/meta/ROD_v1.1_dotbracket.tsv | head | python /watson/rnalab/robedgar/src/py/columns.py
                     sample_name  length     mfe
  DRR000346_804025_circle_804025     240   -73.0
DRR000426_1011226_circle_1011226     262  -161.7
      DRR000601_6173_circle_6173     306   -94.9
  DRR000836_152450_circle_152450     284  -113.6
    DRR000836_16226_circle_16226     349  -146.7
  DRR000836_162517_circle_162517     362  -155.3
    DRR000836_42908_circle_42908     304  -142.7
      DRR000836_5505_circle_5505     400  -175.8
    DRR000836_66877_circle_66877     636  -306.0

# head ~/rce_assiniboine/big_output/uniques.compos.tsv | head | /watson/rnalab/robedgar/src/py/columns.py
#                                       L   NA   NC   NG  NT
  ROD|DRR235453_664787_circle_664787  327   60  110  114  43
    ROD|DRR213772_88310_circle_88310  327   38  125  105  59
ROD|DRR089976_1105891_circle_1105891  415  120  129  122  44
  ROD|DRR235446_961449_circle_961449  363   60  140  112  51
 ROD|ERR1532545_262059_circle_262059  239   57   68   68  46
  ROD|DRR438535_292438_circle_292438  303   63  119   96  25
  ROD|DRR213769_724440_circle_724440  309   80   82  105  42
    ROD|ERR10083444_3366_circle_3366  273   61   81   88  43
  ROD|DRR274963_285854_circle_285854  327   43  114  110  60
  ROD|DRR299359_553199_circle_553199  252   50   73   65  64
'''

home = "/watson/rnalab/robedgar"
composfn = home + "/rce_assiniboine/big_output/uniques.compos.tsv"
mfefn = home + "/assiniboine/meta/ROD_v1.1_dotbracket.tsv"

label2compos = {}
for line in open(composfn):
	line = line.rstrip().replace("ROD|", "")
	flds = line.split('\t')
	label = flds[0]
	label2compos[label] = '\t'.join(flds[1:])

found = 0
notfound = 0
f = open(mfefn)
f.readline()
for line in f:
	line = line.rstrip()
	flds = line.split('\t')
	label = flds[0]
	compos = label2compos.get(label)
	if not compos:
		notfound += 1
		continue
	found += 1
	s = '\t'.join(flds[0:3])
	sys.stdout.write(s + '\t' + compos + "\n")
sys.stderr.write("Found %d, not found %d\n" % (found, notfound))
