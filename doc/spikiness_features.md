# Spikiness features

A rod is one path of helices between two hairpin loops. A spiky drawing has
either pairs that cannot lie on any such path, or several composite stems, or
a junction whose spokes are similar in length. The three percentages below
separate those cases. Each is a percent of paired bases, from 0 to 100, and
each is 0 when the circle has no pairs. They are identity features: they are
not passed through `log1p`.

`_rod_percentages` in `py/dotbracket_circular_features.py` builds the rod
graph. Nodes are hairpin and bifurcation loops. An edge is the helix stretch
that joins two of those loops through any number of bulges and interior
loops. `pRod` is the paired-base percent of the hairpin-to-hairpin path with
the most pairs. On an unbranched fold that path contains every pair, so
`pRod`, `pS`, and `pLongestRod` are equal and `nBranch` is 0.

## pOffRod

```
pOffRod = 100 * (nr_paired_bases - pairs_on_best_path) / nr_paired_bases
```

This is the share of pairs the best rod leaves behind. A bumpy rod scores 0,
because every helix lies on one path even when interior loops break it into
several composite stems. A three-arm star scores about 33, because the best
path takes two arms. A short twig on a long path scores low, and the path is
allowed to bend through junctions.

## pSplitStem

```
pSplitStem = 100 * (1 - longest_composite_stem_nr_paired_bases / nr_paired_bases)
```

Composite stems bridge gaps of at most 5 unpaired bases on a strand, the same
`max_unpaired_gap` as the other stem features. This percentage ignores
branching. Two equal composite stems score 50, whether they form a rod or a
star. One composite stem, however bulged, scores 0. A rod whose interior loop
is longer than 5 therefore scores high here and 0 on `pOffRod`.

## pJunction

At each bifurcation, take the paired length of every incident rod edge. Drop
the heaviest spoke and divide the rest by the sum of the spokes. Average
those local scores, weighting each junction by the paired bases that meet
there. No bifurcation scores 0. A junction with spoke lengths 100, 1, 1
scores about 2. A symmetric three-way junction scores about 67. An edge
between two bifurcations is counted at both ends: this is a property of the
junctions, not a partition of the pairs.

## Complements

The discriminator weights are nonnegative, so a feature that is higher in decoys cannot receive a negative weight. The complement is that negative weight written as a positive one. Each complement is 0 when nothing is paired, not 100, because the source value 0 is the unpaired convention rather than a measured rod.

```
pOnRod = 100 - pOffRod
pLongestPairShare = 100 - pSplitStem
pDominantSpoke = 100 - pJunction
```

`pOnRod` is the share of pairs on the best path. `pLongestPairShare` is the share inside the longest composite stem. `pDominantSpoke` is the share of junction spokes kept by the heaviest spoke. A hairpin scores 100 on all three. The three-arm star scores 66.7, 33.3, and 33.3.

## Where the three disagree

| fold | `pOffRod` | `pSplitStem` | `pJunction` |
|---|---:|---:|---:|
| one hairpin, `(((...)))` | 0 | 0 | 0 |
| two equal stems separated by six dots, `((...))......((...))` | 0 | 50 | 0 |
| three equal arms, `(((...)))(((...)))(((...)))` | 33.3 | 66.7 | 66.7 |

`_check_rod_features` locks these values.
