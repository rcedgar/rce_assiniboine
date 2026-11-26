import sys
from dotbracket_circular_features import *

sample_ss = "((..(((.(((...))))))..))"
feats = extract_circular_rna_features(
	sample_ss, max_unpaired_gap=5, train_median_bulge_len=1.0
)
feature_names = sorted(list(feats.keys()))

# header
print('#label\t' + '\t'.join(feature_names))

for line in open(sys.argv[1]):
	flds = line[:-1].split('\t')
	label = flds[0]
	seq = flds[1]
	ss = flds[2]
	feats = extract_circular_rna_features(
		ss, max_unpaired_gap=5, train_median_bulge_len=1.0
	)
	s = label
	for feature_name in feature_names:
		s += "\t%.3g" % feats[feature_name]
	print(s)