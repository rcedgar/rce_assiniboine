# ROD features vs RCE features

Comparison of the features named in [ROD_features.md](ROD_features.md) with the features returned by `extract_circular_rna_features` in [py/dotbracket_circular_features.py](../py/dotbracket_circular_features.py).

ROD lists 22 features: 6 of sequence composition, 9 of structure composition, and 7 of rod-likeness. RCE returns 24 floats. Two ROD features equal two RCE features. The others are a shared idea with different rules, or they exist on only one side.

RCE definitions follow the function body. The module comment is shorter than the code on the bulge features; those two points are under [Bulges and loop bases](#bulges-and-loop-bases).

The longer ROD inventory in `notes/rod_screen_feature_description.md` (`length`, `mfe`, ORFs, library metadata, enrichment scores) is outside this comparison.

## What each side is looking at

- ROD sequence features are computed from the nucleotide string. RCE reads a dot-bracket of `(`, `)`, and `.` only.
- ROD structure composition is a linear reading of that string. A hairpin is one top-level opening and closing of brackets, and nested branches belong to it. Every base is paired (`n_bound`), unpaired inside a hairpin (`n_loop`), or unpaired outside every pair (`n_unbound`), and those three counts sum to the contig length.
- ROD rod-likeness reads the same string as a circle: a tree of loops joined by base pairs. The bases outside every pair count as one more loop, so two hairpins side by side are one stem that runs across position 1. A loop with three or more base pairs at its edge is a branch point. An unbranched stem runs from one end loop or branch point to the next, and the pairs along it are joined only by bulges and interior loops.
- RCE stem and bulge features are circular from the start. Index steps wrap with modulo the length, so a stem can cross the cut between the last base and the first. A stem stops when the next paired base on one strand is not the partner of the next paired base on the other strand, which is a multiloop.

## Equivalent

The same number on every valid dot-bracket. `length` in the ROD formulas is `nr_bases`.

| ROD | RCE | Quantity |
|:----|:----|:---------|
| `n_bound` | `nr_paired_bases` | Number of `(` and `)`. |
| `pct_bound` | `pct_paired_bases` | `100 * paired / length`. |

## Related, not equivalent

### Unbranched stem and longest composite stem

| ROD | RCE |
|:----|:----|
| `unbranched_length` | `longest_composite_stem_total_nr_bases` |
| `pct_unbranched` | `longest_composite_stem_pct_total_bases` |
| `unbranched_bound` | `longest_composite_stem_nr_paired_bases` |
| `pct_unbranched_bound` | `longest_composite_stem_pct_paired_bases` |

Both describe a run of stacked pairs interrupted by bulges and interior loops, and both can continue that run across position 1. The values differ for four reasons.

**Gap.** ROD keeps joining across a bulge or interior loop of any size, and stops at an end loop or a branch point. RCE keeps joining only while each strand has at most `max_unpaired_gap` unpaired bases between the pairs. The default is 5. A longer interruption ends the RCE stem, and those unpaired bases are absent from the composite-stem length and from the bulge totals.

**End loops.** `unbranched_length` counts the paired bases, the loops between them, and both end loops (a hairpin loop, the exterior loop, or a branch-point loop at the end of the run). `longest_composite_stem_total_nr_bases` counts the paired bases plus the unpaired bases that sit between successive pairs. The apical hairpin loop and the junction or exterior dots are left out. In the code those between-pair bases are `loop_segments`, filled only by the inward step from one pair to the next.

**Which stem wins.** ROD measures `unbranched_length` on the stem with the most bases, and `unbranched_bound` is the paired-base count of that same stem. RCE ranks stems by the number of base pairs (`len(pairs)`) and reports that stem. A stem that is longer because of its loops can win for ROD while a shorter, more fully paired stem wins for RCE.

**Denominator.**

- `pct_unbranched` and `longest_composite_stem_pct_total_bases` are both `100 * stem bases / contig length`, on the two base sets above.
- `pct_unbranched_bound` is `100 * paired bases in that stem / contig length`.
- `longest_composite_stem_pct_paired_bases` is `100 * paired bases / (paired + between-pair unpaired)` inside the winning stem.

### Unpaired bases

`nr_unpaired_bases` equals `n_loop + n_unbound`: every `.`. ROD stores the two parts in separate columns.

- `n_loop` is unpaired bases enclosed by at least one pair.
- `n_unbound` is unpaired bases before, after, or between hairpins in the linear reading.

A base that the linear reading calls unbound can lie on an RCE stem that crosses the origin. That base is then a bulge or interior-loop base in the RCE features, and an exterior base in `n_unbound`.

### Longest hairpin and the stems

`longest_hairpin` is the length of the longest top-level hairpin, from its opening bracket to its closing bracket. The span includes nested branches and the apical loop. `pct_longest_hairpin` is that length as a percentage of the contig.

Two hairpins placed side by side are two hairpins in this count. In the circular reading they are one stem that crosses position 1. The hairpin span therefore differs from `unbranched_length`, from `longest_composite_stem_total_nr_bases`, and from `longest_stem_total_nr_bases`.

### Strict stem

The four `longest_stem_*` features are a helix with gap 0: stacked pairs with no unpaired base between them. ROD has no feature of that helix.

The stem contains only paired bases, so `longest_stem_total_nr_bases` equals `longest_stem_nr_paired_bases`. `longest_stem_pct_paired_bases` is 100 when the structure has a pair, and 0 when it has none. `longest_stem_pct_total_bases` is that perfect helix as a percentage of the contig length. A gap of 0 is the special case of a composite stem, so `longest_stem_nr_paired_bases` is at most `longest_composite_stem_nr_paired_bases`. `pct_unbranched_bound` divides the pairs of a different stem by the contig length.

### Bulges and loop bases

RCE bulge features count per-strand runs of dots that sit between the pairs of a composite stem (each strand at most `max_unpaired_gap`). A run on one strand only is a true bulge (`nr_bases_true_bulges`). Runs on both strands of the same step are a symmetric interior loop (`nr_bases_symmetric_loops`).

`n_loop` counts every unpaired base inside at least one pair: hairpin loops, multiloop interiors, interior loops longer than the gap, and the short bulges. Hairpin-loop dots and exterior dots stay out of the RCE bulge totals. The module comment describes `max_bulge_len` as the longest run of consecutive dots in the string. The code takes the longest of those per-strand composite-stem runs.

`pct_loop` is `n_loop` over the contig length. `pct_bulge`, `pct_true_bulges`, and `pct_symmetric_loops` are percentages of smaller sets. `pct_bulge` splits into the true-bulge and symmetric-loop percentages.

## ROD only

### Sequence composition

These need the nucleotide string. RCE has no sequence input.

| Feature | What it measures |
|:--------|:-----------------|
| `gc_content` | G + C bases as a percentage of the contig length. |
| `keto_skew` | Absolute keto skew, `abs(GT - AC) / (GT + AC)`. |
| `entropy` | Shannon entropy of the base composition of the whole contig (bits). |
| `minEntropy20` | Lowest Shannon entropy of any 20 nt window. Windows do not wrap. |
| `meanEntropy20` | Mean Shannon entropy of all 20 nt windows. |
| `deltaEntropy20` | `minEntropy20 - entropy`. |

### Linear structure

RCE has no column for these. `n_loop` and `n_unbound` are the two parts of `nr_unpaired_bases`, described under [Unpaired bases](#unpaired-bases). `longest_hairpin` is compared with the stems under [Longest hairpin and the stems](#longest-hairpin-and-the-stems).

| Feature | What it measures |
|:--------|:-----------------|
| `n_loop`, `pct_loop` | Unpaired bases enclosed by at least one pair: count, and percentage of length. |
| `n_unbound`, `pct_unbound` | Unpaired bases outside every pair: count, and percentage of length. |
| `n_hairpins` | Number of top-level hairpins. |
| `longest_hairpin`, `pct_longest_hairpin` | Longest top-level hairpin span: nucleotides, and percentage of length. |

### Rod shape

RCE has no column for the backbone or for branch points. The four `unbranched_*` features are the near-miss in [Unbranched stem and longest composite stem](#unbranched-stem-and-longest-composite-stem), so they are omitted here.

| Feature | What it measures |
|:--------|:-----------------|
| `pct_backbone` | Base pairs on the backbone as a percentage of all base pairs. 100 for a rod. |
| `n_branch_points` | Number of loops with three or more pairs at the edge. 0 for a rod. |
| `pct_branch` | Bases in side branches as a percentage of the contig length. |

## RCE only

### Length and unpaired total

| Feature | What it measures |
|:--------|:-----------------|
| `nr_bases` | Length of the dot-bracket. ROD percentages use this length and [ROD_features.md](ROD_features.md) does not name it. |
| `nr_unpaired_bases` | Every `.`. Equal to `n_loop + n_unbound`, which ROD stores separately. |

### Strict stem

Gap 0. See [Strict stem](#strict-stem).

| Feature | What it measures |
|:--------|:-----------------|
| `longest_stem_total_nr_bases` | Bases in the longest perfect helix. Equals the paired count, because the gap is 0. |
| `longest_stem_nr_paired_bases` | Paired bases in that helix. |
| `longest_stem_pct_paired_bases` | `100 * paired / total` inside the helix: 100 if any pair exists, otherwise 0. |
| `longest_stem_pct_total_bases` | `100 * helix bases / nr_bases`. |

The winning helix is the one with the most pairs. Ties keep the first stem in 5' order of its pairs.

### Longest composite stem

Near-miss with `unbranched_*`. See [Unbranched stem and longest composite stem](#unbranched-stem-and-longest-composite-stem). The winner is the stem with the most pairs. A tie keeps the first stem in 5' order, as for the strict stem.

| Feature | What it measures |
|:--------|:-----------------|
| `longest_composite_stem_total_nr_bases` | Paired bases plus between-pair unpaired bases of the stem with the most pairs. |
| `longest_composite_stem_nr_paired_bases` | Paired bases in that stem. |
| `longest_composite_stem_pct_paired_bases` | Those paired bases as a percentage of the stem's own bases. |
| `longest_composite_stem_pct_total_bases` | The stem's bases as a percentage of `nr_bases`. |

### Bulges

Twelve features, all taken from the per-strand dot runs inside composite stems. See [Bulges and loop bases](#bulges-and-loop-bases).

Two details are in the code and absent from the module comment:

- An excess bulge is a run longer than `train_median_bulge_len`. The default is 1, so a run of 2 or more is excess at the default.
- `bulge_to_stem_ratio` is `total_nr_bulge_bases / nr_paired_bases` (0 when there are no pairs). The comment stops at `total_nr_bulge_bases/`.

| Feature | What it measures |
|:--------|:-----------------|
| `max_bulge_len` | Longest per-strand run inside a composite stem. |
| `nr_excess_bulges` | How many of those runs exceed `train_median_bulge_len`. |
| `total_nr_bulge_bases` | Bases in all such runs. |
| `pct_bulge` | `100 * total_nr_bulge_bases / nr_bases`. |
| `excess_bulge_nr_bases` | Bases that sit in excess runs. |
| `pct_excess_bulges` | `100 * excess_bulge_nr_bases / nr_bases`. |
| `bulge_to_stem_ratio` | `total_nr_bulge_bases / nr_paired_bases`. |
| `bulge_length_L2_norm` | `sqrt(sum of squared run lengths)`. |
| `nr_bases_true_bulges` | Bases in one-sided runs. |
| `nr_bases_symmetric_loops` | Bases in two-sided runs. |
| `pct_true_bulges` | `100 * nr_bases_true_bulges / nr_bases`. |
| `pct_symmetric_loops` | `100 * nr_bases_symmetric_loops / nr_bases`. |

A structure with no composite-stem loop segment leaves every bulge feature at 0. `total_nr_bulge_bases` equals `nr_bases_true_bulges + nr_bases_symmetric_loops`.

## Counts

| Class | ROD | RCE |
|:------|----:|----:|
| Equivalent | 2 | 2 |
| Unbranched stem / composite stem | 4 | 4 |
| Linear structure with no RCE column | 7 | 0 |
| Sequence composition | 6 | 0 |
| Backbone and branch points | 3 | 0 |
| Length, unpaired total, strict stem, bulges | 0 | 18 |
| **Total** | **22** | **24** |

The 7 linear-structure features are `n_loop`, `pct_loop`, `n_unbound`, `pct_unbound`, `n_hairpins`, `longest_hairpin`, and `pct_longest_hairpin`. The 18 RCE-only features are `nr_bases`, `nr_unpaired_bases`, the four `longest_stem_*` features, and the twelve bulge features. The four `longest_composite_stem_*` features are the composite-stem row, not that 18.
