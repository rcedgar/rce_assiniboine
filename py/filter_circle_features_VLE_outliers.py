#!/usr/bin/env python3

import sys

infn = "../big_output/controls_circle_features.tsv"
outfn = "../big_output/controls_circle_features_filtered.tsv"

f = open(infn)
outf = open(outfn, "w")
raw_header = f.readline()
outf.write(raw_header)
hdr = raw_header[1:-1]
flds = hdr.split('\t')
name2col = {}
for i in range(len(flds)):
    fld = flds[i]
    name2col[fld] = i

pct_bulge_col = name2col["pct_bulge"]
pct_symmetric_loops_col = name2col["pct_symmetric_loops"]
pct_excess_bulges_col = name2col["pct_excess_bulges"]
bulge_to_stem_ratio_col = name2col["bulge_to_stem_ratio"]
bulge_length_L2_norm_col = name2col["bulge_length_L2_norm"]

npct_bulge = 0
npct_symmetric_loops = 0
npct_excess_bulges = 0
nbulge_to_stem_ratio = 0
nbulge_to_stem_L2_norm = 0
nfiltered = 0
navsunviroid = 0

for line in (f):
    flds = line[:-1].split('\t')
    label = flds[0]
    parts = label.split("_", 2)
    if len(parts) >= 2 and parts[1] == "avsunviroid":
        navsunviroid += 1
        continue
    pct_bulge = float(flds[pct_bulge_col])
    pct_symmetric_loops = float(flds[pct_symmetric_loops_col])
    pct_excess_bulges = float(flds[pct_excess_bulges_col])
    bulge_to_stem_ratio = float(flds[bulge_to_stem_ratio_col])
    bulge_to_stem_L2_norm = float(flds[bulge_length_L2_norm_col])
    if label.startswith("VLE_"):
        if pct_bulge < 12:
            npct_bulge += 1
            continue
        if pct_symmetric_loops < 10:
            npct_symmetric_loops += 1
            continue
        if pct_excess_bulges < 8:
            npct_excess_bulges += 1
            continue
        if bulge_to_stem_ratio < 0.175:
            nbulge_to_stem_ratio += 1
            continue
        if bulge_to_stem_L2_norm < 10:
            nbulge_to_stem_L2_norm += 1
            continue
        nfiltered += 1
    outf.write(line)
outf.close()

sys.stderr.write(f"{npct_bulge=}\n")
sys.stderr.write(f"{npct_symmetric_loops=}\n")
sys.stderr.write(f"{npct_excess_bulges=}\n")
sys.stderr.write(f"{nbulge_to_stem_ratio=}\n")
sys.stderr.write(f"{nbulge_to_stem_L2_norm=}\n")
sys.stderr.write(f"{navsunviroid=}\n")
sys.stderr.write(f"{nfiltered=}")
