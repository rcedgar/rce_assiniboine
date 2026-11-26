# ROD v2 features vs RCE features

Comparison of the features named in [ROD_features_v2.md](ROD_features_v2.md) with the 24 floats returned by `extract_circular_rna_features` in [py/dotbracket_circular_features.py](../py/dotbracket_circular_features.py).

RCE reads one dot-bracket of `(`, `)`, and `.`. ROD v2 is the dictionary of the contig tables, the reverse-complement fold, and the ORF, library, score, UMAP, and cluster tables. Definitions of the RCE features follow the function body. The module comment is shorter than the code on the bulge features; those points are under [Bulges and interior loops](#bulges-and-interior-loops).

The earlier comparison against [ROD_features.md](ROD_features.md) is [ROD_features_vs_rce_features.md](ROD_features_vs_rce_features.md). What that note treated as current columns, and v2 has since dropped, is listed at the end.

## What each side is looking at

- ROD sequence features use the nucleotide string. The fold features use that string together with the dot-bracket. RCE never sees the bases.
- Both readings of a fold are circular. ROD: the bases outside every pair are one more loop, two hairpins side by side are one stem across position 1, and rotating the start leaves the features unchanged. RCE: index steps wrap with modulo the length, so a stem can cross the cut between the last base and the first.
- ROD gives every base one of five states. A paired base is a stack (S). An unpaired base takes the state of its loop, classed by how many pairs meet that loop: a hairpin loop (H, one pair), an interior loop (I, two pairs, unpaired on both sides), a bulge (B, two pairs, unpaired on one side), or a bifurcation loop (F, three or more pairs). `pS + pH + pI + pB + pF` is 100 when the structure has a pair. All five are 0 when it has none.
- RCE stops a stem when the next paired base on one strand is not the partner of the next paired base on the other, and it joins two helices only when each strand has at most `max_unpaired_gap` unpaired bases between them. The default gap is 5. Hairpin-loop bases and bifurcation-loop bases stay out of the stem totals and out of the bulge totals.

Every structure feature below is also stored for the reverse complement, under the same name with `_rc`, from `dot_bracket_rc`. RCE has no second set of outputs. The same 24 numbers come back if that string is passed in instead.

## Equivalent

The same number on every contig whose dot-bracket has one character per base, which is how `dot_bracket` is stored.

| ROD | RCE | Quantity |
|:----|:----|:---------|
| `length` | `nr_bases` | Number of positions. ROD counts the FASTA sequence. RCE counts the dot-bracket. |
| `pS` | `pct_paired_bases` | Paired bases, `(` or `)`, as a percentage of that length. `pS` is the old `pct_bound`. Both are 0 when the structure has no pair. |

Two further RCE columns are fixed by those two. ROD does not store them.

- `nr_paired_bases = length * pS / 100`
- `nr_unpaired_bases = length * (100 - pS) / 100`

The second identity uses `100 - pS`, which is the unpaired percentage both when the five states sum to 100 and when a pair-free structure has set `pS` to 0. The sum `pH + pI + pB + pF` equals that unpaired percentage only in the first case. On a pair-free structure the four loop percentages are 0 and `nr_unpaired_bases` is `length`.

`pS_rc` is `pct_paired_bases` of the reverse-complement fold, not a second RCE feature.

## Related, not equivalent

### Bulges and interior loops

| ROD | RCE | Shared idea |
|:----|:----|:------------|
| `pB` | `pct_true_bulges` | Unpaired bases on one side of a two-pair loop, as a percentage of length. |
| `pI` | `pct_symmetric_loops` | Unpaired bases on both sides of a two-pair loop, as a percentage of length. |
| `pB + pI` | `pct_bulge` | Those two sets together, as a percentage of length. |

The base sets match when every two-pair loop has at most `max_unpaired_gap` unpaired bases on each side. A longer loop stays in `pB` or `pI` and drops out of the RCE bulge features, because the inward step refuses the bridge and the stem ends. The default gap is 5.

RCE records each side of a bridged step as its own run. One positive side is a true bulge. Two positive sides are a symmetric loop. `total_nr_bulge_bases` equals `nr_bases_true_bulges + nr_bases_symmetric_loops`, and `pct_bulge` is that sum over `nr_bases`.

`pH` and `pF` have no RCE column. Hairpin-loop bases and bifurcation-loop bases are outside the RCE bulge totals. The module comment describes `max_bulge_len` as the longest run of consecutive dots in the string. The code takes the longest per-strand run inside a composite stem.

The other bulge columns are RCE only. Two of them are specified in the code and not in the comment: an excess run is longer than `train_median_bulge_len` (default 1, so a run of 2 or more is excess), and `bulge_to_stem_ratio` is `total_nr_bulge_bases / nr_paired_bases` (0 when there are no pairs).

| RCE feature | What it measures |
|:------------|:-----------------|
| `max_bulge_len` | Longest per-strand run inside a composite stem. |
| `nr_excess_bulges` | How many of those runs exceed `train_median_bulge_len`. |
| `total_nr_bulge_bases` | Bases in all such runs. |
| `excess_bulge_nr_bases` | Bases that sit in excess runs. |
| `pct_excess_bulges` | `100 * excess_bulge_nr_bases / nr_bases`. |
| `bulge_to_stem_ratio` | `total_nr_bulge_bases / nr_paired_bases`. |
| `bulge_length_L2_norm` | `sqrt(sum of squared run lengths)`. |
| `nr_bases_true_bulges` | Bases in one-sided runs. The percentage form is the `pB` partner above. |
| `nr_bases_symmetric_loops` | Bases in two-sided runs. The percentage form is the `pI` partner above. |

A structure with no bridged loop segment leaves every bulge feature at 0.

### Backbone, longest stretch, and the RCE stems

ROD's backbone is the path from one hairpin loop to another that holds the most base pairs. A tie goes to the path with the longest stretch. Branch points on that path cut it into unbranched stretches: runs of pairs joined by bulges and interior loops. Side branches hang off the path.

| ROD | What it measures | Nearest RCE feature |
|:----|:-----------------|:--------------------|
| `pRod` | Paired bases of the whole backbone, as a percentage of length. Equal to `pS` for a rod. 0 when there is no pair. | None. A backbone can be several stretches. RCE reports one stem. |
| `pLongestRod` | Paired bases of the longest unbranched stretch of that backbone, as a percentage of length. `pLongestRod <= pRod <= pS`, with equality for a rod that has no branch point. | `longest_composite_stem_nr_paired_bases`, divided by `nr_bases`. RCE does not return that percentage. |
| `nBranch` | Number of bifurcation loops. 0 for a rod. | None. RCE uses a bifurcation as a place to stop, and does not count it. |

`pRod` replaces `pct_backbone`. The old column was backbone pairs as a percentage of all pairs; `pRod = pct_backbone * pS / 100`.

`pLongestRod` replaces `pct_unbranched_bound` and is narrower. The old column took the longest unbranched stem anywhere, including a side branch. `pLongestRod` stays on the backbone.

Against the four `longest_composite_stem_*` columns the stretch still differs:

- **Gap.** ROD joins a bulge or interior loop of any size. RCE joins it only while each strand has at most 5 unpaired bases. A longer interruption ends the RCE stem and removes those bases from the bulge totals.
- **Where the stem may lie.** `pLongestRod` is a stretch of the backbone. RCE keeps the stem with the most pairs anywhere, including a side branch. The winner is that stem; a tie keeps the first stem in 5' order.
- **What is counted.** `pLongestRod` is paired bases over the contig length. `longest_composite_stem_nr_paired_bases` is the paired count of the RCE winner. `longest_composite_stem_total_nr_bases` adds the unpaired bases between those pairs and leaves out both hairpin loops and any bifurcation loop. `longest_composite_stem_pct_total_bases` is that total over the contig length. `longest_composite_stem_pct_paired_bases` is the paired bases over the stem's own bases.

On a rod with no branch point and no gap above 5, `pLongestRod` equals `100 * longest_composite_stem_nr_paired_bases / nr_bases`. `longest_composite_stem_pct_total_bases` is larger by the between-pair unpaired bases. A gap above 5, or a side branch with more pairs than any backbone stretch, separates the two.

### Strict stem

The four `longest_stem_*` features are a helix with gap 0. ROD has no feature of that helix. The stem contains only paired bases, so `longest_stem_total_nr_bases` equals `longest_stem_nr_paired_bases`. `longest_stem_pct_paired_bases` is 100 when the structure has a pair, and 0 when it has none. `longest_stem_pct_total_bases` is that helix as a percentage of length.

A gap of 0 is the special case of a composite stem, so `longest_stem_nr_paired_bases` is at most `longest_composite_stem_nr_paired_bases`. The same 5' tie-break applies. `pLongestRod` still includes pairs that a bulge or an interior loop separates, and it ignores a perfect helix that lies on a side branch.

## ROD only

### Loop states, keto skew, and the other strand

| Feature | What it measures |
|:--------|:-----------------|
| `pH` | Hairpin-loop bases, loops with one pair, as a percentage of length. On a circle a single helix has two such loops: the apical loop and the exterior loop across position 1. |
| `pF` | Bifurcation-loop bases, loops with three or more pairs, as a percentage of length. |
| `bKeto` | Absolute keto skew of the paired bases, `abs(GT - AC) / (GT + AC)`, from 0 to 1. A G-C or A-U pair contributes one keto base and a G-U pair two, so this is the share of G-U pairs among the pairs. 0 when there is no pair. |
| `uKeto` | Absolute keto skew of the unpaired bases, from 0 to 1. |
| `pS_rc`, `pH_rc`, `pI_rc`, `pB_rc`, `pF_rc`, `pRod_rc`, `nBranch_rc`, `pLongestRod_rc`, `bKeto_rc`, `uKeto_rc` | The structure block, including the partners of RCE features, computed on `dot_bracket_rc`. |

`bKeto` and `uKeto` replace `keto_skew` of the whole contig. The LSC table still keeps `keto_skew`; the fold tables do not. Both need the sequence and the pairing. RCE has neither.

### Sequence, graph, and fold energy

RCE does not take a sequence, a header, or an energy.

| Feature | What it measures |
|:--------|:-----------------|
| `contig`, `header`, `sequence` | Name, FASTA header, and nucleotide string. `dot_bracket` is computed from `sequence`. |
| `length_mod3` | `TRUE` when `length` is a multiple of 3. |
| `kmer_coverage` | BCALM2 `ka` abundance. |
| `gc_content` | G + C as a percentage of length. |
| `entropy` | Shannon entropy of the whole-contig base composition (bits). |
| `minEntropy20` | Lowest Shannon entropy of any 20 nt window. Windows do not wrap. |
| `meanEntropy20` | Mean Shannon entropy of those windows. |
| `deltaEntropy20` | `minEntropy20 - entropy`. |
| `circular` | Assembly-graph self-closure, from the header edges. |
| `mfe`, `mfent`, `mfentgc`, `mfeZ` | Fold energy, energy per nucleotide, energy relative to GC, and the Z-score of `mfe` against shuffled sequences. `mfentgc` is calculated for figures and is not saved. |
| `mfe_rc`, `mfent_rc` | The same energies for the reverse complement. |
| `dot_bracket`, `dot_bracket_rc` | The structures themselves. Either string can be the RCE input. Neither is an RCE output. |

Each Parquet table also has `row`, the contig's row number.

### Open reading frames, library, scores

| Feature | What it measures |
|:--------|:-----------------|
| `orf_length`, `orf_strand`, `orf_code`, `n_orfs`, `orf_pct` | Longest ORF and how many ORFs the contig has. `orf_pct` can pass 100 because ORFs are called on the circle. |
| `acc`, `organism`, `assay_type`, `librarysource`, `biosample`, `bioproject` | SRA library of the contig. |
| `group`, `z`, `z_norm`, `fc`, `fc_norm`, `sig` | RNA-Seq enrichment scores summed over 14 binned features. |

### Maps and clusters

These are coordinates, labels, or summaries of many contigs. None is a function of one dot-bracket.

Per contig: `umap_x`, `umap_y`, `umap_cluster`, `rna_umap_x`, `rna_umap_y`, `rna_umap_cluster`, `RNA_cluster`, `orf_first`, `orf_cluster70`, `orf_cluster30`, `orf_cluster_members`, `orf_cluster_fold`, `orf_cluster_enriched`, `orf_model`.

Per ORF (`orf.df`): `category`, `contig_row`, `orf_nr`, `orf_start`, `orf_end`, `strand`, `frame`, `orf_code`, `aa_length`, `cluster70`, `cluster30`, `is_centroid`, `identity`, `member_cover`, `centroid_cover`.

Per ORF cluster (`orf.cluster70`, `orf.cluster30`): `centroid`, `centroid_category`, `centroid_contig_row`, `centroid_orf_nr`, `centroid_aa_length`, `primary_orf`, `cluster30`, `clusters70`, `members`, `n_rod`, `n_decoy`, `n_vle`, `contigs`, `n_test`, `n_reference`, `expected`, `fold`, `z`, `fold_circular`, `z_circular`, `enriched`, `meanLength`, `meanGC`, `meanMFE`, `sdLength`, `meanRod`, `sdRod`, `n_libraries`, `n_bioprojects`, `percent_RNA`, `organism`, `organism_pct`, `organisms`, `umap_cluster`, `umap_cluster_pct`, `umap_enriched`, `umap_enriched_pct`, `umap_enriched_fold`, `rna_umap_cluster`, `rna_umap_cluster_pct`, `model`, `model_evalue`, `model_bitscore`, `model_clusters70`. `meanRod` and `sdRod` summarize `pRod`.

Per UMAP cluster (`umap_clusters`): `contigs`, `x`, `y`, `pct_circular`, `n_test`, `n_reference`, `expected`, `fold`, `z`, `fold_circular`, `z_circular`, `enriched`, `n_vle`, `fold_vle`, `z_vle`, `n_decoy`, `fold_decoy`, `z_decoy`.

Per RNA-cluster (`umap_rna_clusters`): `contigs`, `x`, `y`, `n_test`, `n_reference`, `fold`, `z`, `enriched`, `fold_within`, `n_runs`, `top_run_pct`, `top_project_pct`, `origin`, `origin_pct`, `n_vle`, `n_decoy`, `element`, `element_n`, `element_fold`, `elements`, `RNA_cluster`.

`umap3d` places a sample in three dimensions: `source`, `source_row`, `x`, `y`, `z`.

## RCE only

Stored on the RCE side alone. The circle-level counts are the identities under [Equivalent](#equivalent). The bulge percentages and the composite stem are the near-misses under [Related, not equivalent](#related-not-equivalent).

| Feature | What it measures |
|:--------|:-----------------|
| `nr_paired_bases` | Count of `(` and `)`. Equal to `length * pS / 100`. |
| `nr_unpaired_bases` | Count of `.`. Equal to `length * (100 - pS) / 100`. |
| `longest_stem_total_nr_bases` | Bases in the longest perfect helix. Equals the paired count. |
| `longest_stem_nr_paired_bases` | Paired bases in that helix. |
| `longest_stem_pct_paired_bases` | 100 if any pair exists, otherwise 0. |
| `longest_stem_pct_total_bases` | `100 * helix bases / nr_bases`. |
| `longest_composite_stem_total_nr_bases` | Pairs plus between-pair unpaired bases of the stem with the most pairs. |
| `longest_composite_stem_nr_paired_bases` | Paired bases in that stem. |
| `longest_composite_stem_pct_paired_bases` | Those pairs as a percentage of the stem's own bases. |
| `longest_composite_stem_pct_total_bases` | The stem's bases as a percentage of `nr_bases`. |
| `max_bulge_len`, `nr_excess_bulges`, `total_nr_bulge_bases`, `excess_bulge_nr_bases`, `pct_excess_bulges`, `bulge_to_stem_ratio`, `bulge_length_L2_norm`, `nr_bases_true_bulges`, `nr_bases_symmetric_loops` | Bulge size, excess over the training median, and the counts behind `pct_bulge`, `pct_true_bulges`, and `pct_symmetric_loops`. |

## Counts

Forward-strand structure features, plus `length`, against the 24 RCE outputs.

| Class | ROD | RCE |
|:------|----:|----:|
| Equivalent columns (`length`, `pS`) | 2 | 2 |
| Loop or stem quantities that share an idea and differ (`pB`, `pI`, `pRod`, `pLongestRod`) | 4 | 7 |
| ROD structure columns with no RCE partner (`pH`, `pF`, `nBranch`, `bKeto`, `uKeto`) | 5 | 0 |
| RCE columns with no ROD partner (strict stem, bulge counts and norms, and the two circle counts ROD does not store) | 0 | 15 |
| **Total** | **11** | **24** |

The 7 RCE columns in the shared row are `pct_true_bulges`, `pct_symmetric_loops`, and `pct_bulge` (the partners of `pB` and `pI`) and the four `longest_composite_stem_*` features (the partners of `pLongestRod`). `pRod` is in that row on the ROD side and has no RCE column. The 15 are `nr_paired_bases`, `nr_unpaired_bases`, the four `longest_stem_*` features, and the nine bulge columns in the RCE-only table. `nr_paired_bases` and `nr_unpaired_bases` are listed there because ROD does not store them; their values follow from `length` and `pS`.

The same 10 structure features exist again as `*_rc`. Everything else named in [ROD_features_v2.md](ROD_features_v2.md) is ROD only.

## Dropped since the v1 comparison

v2 records these as removed. They are the columns the previous note matched or contrasted, and they are not features of the current tables.

| Removed column | Where it went |
|:---------------|:--------------|
| `pct_bound` | Renamed `pS`. Still equal to `pct_paired_bases`. |
| `n_bound` | Not stored. The count is `length * pS / 100`, which is `nr_paired_bases`. |
| `n_loop`, `n_unbound`, `pct_unbound` | Removed. A circle has no outside. Unpaired bases are `pH`, `pI`, `pB`, and `pF`. |
| `pct_loop` | Was every unpaired base inside a hairpin. `pH` is hairpin-loop bases only. |
| `n_hairpins`, `longest_hairpin`, `pct_longest_hairpin` | Removed. |
| `pct_backbone` | Replaced by `pRod`, which divides backbone pairs by the contig length. |
| `n_branch_points` | Renamed `nBranch`. |
| `pct_unbranched_bound` | Replaced by `pLongestRod`, backbone only. |
| `unbranched_length`, `pct_unbranched`, `unbranched_bound`, `pct_branch` | Removed. |
| `keto_skew` | Replaced on the fold tables by `bKeto` and `uKeto`. |
