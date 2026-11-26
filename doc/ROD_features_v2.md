# Features of ROD v1.1

Every feature of the 11,308,939 contigs of ROD v1.1 used in this project:
the values read from the original data files and the values calculated from
them by the `Sunburst` package and the notebooks of this repository.

- **Origin** says which is which: `file` for a value read unchanged from a
  data file, `derived` for a value calculated in this project.
- **One row per contig.** Every table has the same rows in the same order
  (the order of `ROD_v1.1.fa`), so row `i` is the same contig in each table.
  The contig names are in `rod.df$contig`.
- **Table** is the Parquet file that holds the feature, in the folder
  `tables/` (`rod.df` is `tables/rod.df.parquet`, and so on). Each file also
  has the row number in a first column, `row`. Several features are also
  copied into `rod.df`.
- **Percentages** run from 0 to 100 and are of the contig length unless
  stated otherwise.

| Data file | Holds |
|:----------|:------|
| `ROD_v1.1.fa` | Contig sequences, with the k-mer coverage and assembly-graph edges in the header. |
| `ROD_v1.1.mfeZ.tsv`, `CONTROL_v260921.mfeZ.tsv` | Z-score of the minimum free energy of each sequence against shuffled sequences. |
| `ROD_v1.1.fold.tsv.zst` | RNAfold prediction of each contig and of its reverse complement: minimum free energy and structure of both. (Until 2026-11-21: `ROD_v1.1_dotbracket.tsv`, with the forward strand only.) |
| `ROD_v1.1.reorf_n60.fa` | Open reading frames called on the contigs, one record per ORF. |
| `ROD_v1.1_SRARunInfo.csv` | Metadata of the SRA run each contig was assembled from. |
| `cluster/` | The ORF clusters of the contigs and the controls (a copy of `s3://assiniboine/cluster/`, made with DIAMOND by `inst/scripts/orf_clusters/`): every ORF with its cluster at 70% protein identity, the 70% clusters grouped again at 30%, and the matches of the 70% centroids to protein models. |

## Contig and sequence

Read from `ROD_v1.1.fa`. Each record is a BCALM2 unitig with a header such as
`>DRR000601_6173_circle_6173 ka:f:9.342 L:-:6173:- L:+:6173:+`.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `contig` | `rod.df` | file | Contig name, `<SRA run>_<unitig id>_circle_<unitig id>`: the header up to the first space. The unitig id is written twice in every name, so figures, tables of the explorer and addresses show the name without the repeat (`DRR000601_6173` for `DRR000601_6173_circle_6173`): `contig_label()` makes it where it is shown, `contig_full_name()` and `match_contig()` go back. The tables keep the name of the file, and the short name is not stored. |
| `header` | `rod.fa` | file | Full header line without the `>`: contig name, k-mer coverage and assembly-graph edges. |
| `sequence` | `rod.fa` | file | Nucleotide sequence of the contig, as written in the FASTA. |
| `length` | `rod.df` | file | Contig length: number of bases of the sequence (nt). |
| `length_mod3` | `rod.df` | derived | `TRUE` when `length` is a multiple of 3: a circle of such a length can be read round and round in one reading frame. |
| `kmer_coverage` | `rod.df` | file | k-mer coverage (`ka`): the number of the header field `ka:f:`, so `ka:f:9.342` is 9.342. BCALM2 reports it for every unitig: the mean abundance of its k-mers in the reads of the library, i.e. how deeply the contig was sequenced. |

## Sequence composition

Computed from the nucleotide sequence by `read_rod_fasta()`. `U` counts as
`T`, case is ignored, and other letters (e.g. `N`) are left out.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `gc_content` | `rod.df` | derived | G + C bases as a percentage of the contig length. |
| `entropy` | `rod.entropy`, `rod.df` | derived | Shannon entropy of the base composition of the whole contig, `-sum(p * log2(p))` over A, C, G, T (bits). 2 when the four bases are equally common, 0 for a homopolymer. |
| `minEntropy20` | `rod.entropy` | derived | Lowest Shannon entropy of any 20 nt window of the contig (bits): the least complex 20 nt stretch. Windows start at every base and do not wrap around the end of the contig. |
| `meanEntropy20` | `rod.entropy` | derived | Mean Shannon entropy of all 20 nt windows (bits). |
| `deltaEntropy20` | `rod.entropy` | derived | `minEntropy20 - entropy` (bits). 0 or below; strongly negative when a contig that uses all four bases has a low-complexity stretch. |

## Assembly graph

Computed from the BCALM2 edges in the FASTA header. Each field
`L:<strand>:<unitig id>:<strand>` is an edge of the compacted de Bruijn graph;
the edges themselves are kept in `header`.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `circular` | `rod.df` | derived | `TRUE` when the contig closes on itself: its header has both self-edges `L:+:<id>:+` and `L:-:<id>:-` and no reverse-complement self-edge (`L:+:<id>:-` or `L:-:<id>:+`). A contig with all four self-edges joins its own reverse complement and is counted as an assembly artefact (`FALSE`). |

## Predicted fold

The RNAfold prediction of each contig, read from `ROD_v1.1.fold.tsv.zst`,
and the stability measures calculated from it. The file has two folds per
contig: of the sequence as it is in the FASTA (the forward strand) and of
its reverse complement.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `mfe` | `rod.fold`, `rod.df` | file | Minimum free energy (MFE) of the predicted fold (kcal/mol). More negative is more stable. |
| `dot_bracket` | `rod.dotbracket` | file | Predicted secondary structure in dot-bracket notation, one character per base: `(` and `)` for the two bases of a pair, `.` for an unpaired base. |
| `mfent` | `rod.fold` | derived | Minimum free energy per nucleotide, `mfe / length` (kcal/mol/nt). MFE grows with length; this is the length-independent measure of fold stability. |
| `mfentgc` | not saved | derived | `log(-mfent / gc_content)` (natural log): fold stability relative to GC content. Calculated for the figures of `02B_explore_features.Rmd` only. |
| `mfeZ` | `rod.fold`, `rod.df` | file | Z-score of `mfe` (forward strand): how far it lies from the MFE of shuffled sequences of the same composition (a model of 1,000 shuffles), in standard deviations. Well below 0 for a sequence that folds more stably than its composition explains. Read from `ROD_v1.1.mfeZ.tsv`, which holds each sequence once: a contig that is not in it has the value of the contig with the identical sequence. |
| `mfe_rc` | `rod.fold_rc`, `rod.df` | file | Minimum free energy of the predicted fold of the reverse complement (kcal/mol). |
| `dot_bracket_rc` | `rod.dotbracket_rc` | file | Predicted structure of the reverse complement, written along the reverse complement: its position `i` is base `length + 1 - i` of the contig. |
| `mfent_rc` | `rod.fold_rc` | derived | `mfe_rc / length` (kcal/mol/nt). |

**Reverse complement.** Every feature of the next section is also
calculated, in the same way, from `dot_bracket_rc`. These are in the table
`rod.fold_rc` under the same name with `_rc` at the end: `pS_rc`,
`pH_rc`, `pRod_rc`, `nBranch_rc`, `bKeto_rc` and so on (the two keto skews
from the bases of the reverse complement). The
control sequences have the same tables (`control.fold_rc`,
`control.dotbracket_rc`), from `CONTROL_v260921.fold.tsv.zst`.

## Structure

Computed from the dot-bracket structure by `dot_bracket_features()`, in one
pass (the two keto skews from the structure and the sequence). The
structure is read as a circular molecule, as a tree of loops joined by base
pairs.

- Every base has one of five **states**. A paired base is in a stack (S).
  An unpaired base has the state of its loop, which is classed by the
  number of base pairs at its edge: a **hairpin loop** (H, one pair) closes
  a helix; a loop with two pairs is an **interior loop** (I, unpaired bases
  on both of its sides) or a **bulge** (B, on one side only); a loop with
  three or more is a **bifurcation loop** (F), a **branch point**.
- The bases outside every pair count as one more loop. A structure with one
  hairpin is a helix with a hairpin loop on either side (the hairpin loop,
  and the two ends of the string joined across position 1); two hairpins
  side by side are one stem that runs across position 1, and the bases
  between and around them are its interior loop or bulge. Where a circular
  contig starts therefore does not change its features, and there is no
  state for "outside every pair".
- The **backbone** is the path from one hairpin loop to another that holds
  the most base pairs; what hangs off it is a side branch.
- The branch points on the backbone cut it into **unbranched stretches**:
  runs of base pairs joined only by bulges and interior loops. When several
  paths hold the most base pairs, the one with the longest stretch is the
  backbone.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `pS` | `rod.fold` | derived | Stack: paired bases, `(` or `)`, as a percentage of the contig length. (`pct_bound` until snapshot `v261124b`.) |
| `pH` | `rod.fold` | derived | Hairpin loop: unpaired bases of the loops with one base pair at their edge, as a percentage of the contig length. (`pct_loop` until `v261124b`.) |
| `pI` | `rod.fold` | derived | Interior loop: unpaired bases of the loops with two base pairs at their edge and unpaired bases on both sides. |
| `pB` | `rod.fold` | derived | Bulge: unpaired bases of the loops with two base pairs at their edge and unpaired bases on one side. |
| `pF` | `rod.fold` | derived | Bifurcation loop: unpaired bases of the loops with three or more base pairs at their edge. `pS + pH + pI + pB + pF` is 100; all five are 0 for a structure without base pairs. |
| `pRod` | `rod.fold` | derived | Rod backbone: paired bases of the backbone as a percentage of the contig length. Equal to `pS` for a rod; 0 for a structure without base pairs. (Until `v261124b`: `pct_backbone`, the base pairs of the backbone as a percentage of all base pairs; `pRod = pct_backbone * pS / 100`.) |
| `nBranch` | `rod.fold` | derived | Number of branch points. 0 for a rod. (`n_branch_points` until `v261124b`.) |
| `pLongestRod` | `rod.fold` | derived | Paired bases of the longest unbranched stretch of the backbone as a percentage of the contig length. `pLongestRod <= pRod <= pS`, with equality for a rod without a branch point. (Until `v261124b`: `pct_unbranched_bound`, the longest unbranched stem wherever it lay, also in a side branch.) |
| `bKeto` | `rod.fold` | derived | Absolute keto skew of the paired bases, `abs(GT - AC) / (GT + AC)` with `GT` and `AC` the numbers of keto (G, T) and amino (A, C) bases among them, 0 to 1. A G-C or A-U pair has one of each and a G-U pair two keto bases, so this is the share of G-U pairs among the base pairs. 0 for a structure without base pairs. |
| `uKeto` | `rod.fold` | derived | Absolute keto skew of the unpaired bases, 0 to 1. In most contigs they lean to the amino bases: the G-U pairs take keto bases out of them. |

`bKeto` and `uKeto` take the place of `keto_skew` (in `rod.df` until
snapshot `v261124b`), the absolute keto skew of the whole contig, which is
the difference of the two parts in most contigs and is not kept. The LSC
table (`lsc.df`), which has no folds, keeps its `keto_skew`.

Until snapshot `v261122` the tables also held counts of bases (`n_bound`,
`n_loop`, `n_unbound`, `unbranched_length`, `unbranched_bound`), hairpin
measures (`n_hairpins`, `longest_hairpin`, `pct_longest_hairpin`),
`pct_unbound`, `pct_unbranched` and `pct_branch`, and `pct_loop` meant every
unpaired base inside a hairpin. These were dropped: inside and outside a
hairpin are not different places on a circle.

## Open reading frames

Read from `ROD_v1.1.reorf_n60.fa`, which has one record per ORF of 180 nt (60
codons) or longer, with a header such as
`>DRR000346_804025_circle_804025 ka:f:23.410 L:-:804025:- L:+:804025:+ orf1 code1 240-479/240(-1)Cem`:
the header of the contig as in the FASTA, ORF number, genetic code table,
`start-end`, contig length, `(strand, frame)` and an optional ORF code. The
protein sequences are not kept. Since snapshot `v261124` the file is the one
of 2026-11-23, which has the records of every contig (the file before it
held each sequence once, and a contig with the sequence of another had no
record).

ORFs are called on the contig as a circle. A reading frame without a stop
codon has an ORF code and is reported whole: once around the circle when the
contig length is divisible by 3 (`Ce` with an ATG, `Cem` without one), three
times around when it is not (`Co`, `Com`). The control ORF file follows the
same rule, and `find_orfs()` of the package reproduces both files. (Up to
snapshot `v261123` such a frame of a ROD contig was counted from its first
ATG to the end of the contig only.)

The tables hold the **longest ORF** of each contig: the record with the
largest length, and on ties the first in the file. A contig without a record
has no ORF of 180 nt or longer.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `orf_length` | `rod.orf`, `rod.df` | file | Length of the longest ORF, `end - start + 1` (nt); 0 if the contig has no ORF. |
| `orf_strand` | `rod.orf`, `rod.df` | file | Strand of the longest ORF, as a number: 1 for `+`, -1 for `-`, 0 for a contig without an ORF. |
| `orf_code` | `rod.orf` | file | ORF code of the longest ORF (`Cem`, `Ce`, `Co`, `Com`), `.` when its header has none; `NA` if the contig has no ORF. |
| `n_orfs` | `rod.orf` | derived | Number of ORFs called on the contig; 0 if none. |
| `orf_pct` | `rod.orf` | derived | Longest ORF as a percentage of the contig length. ORFs are called on the contig as a circle, so the value passes 100 when an ORF wraps around it. |

## Library metadata

Read from `ROD_v1.1_SRARunInfo.csv`, which has one row per SRA run. Each
contig takes the values of the run it was assembled from; they are `NA` when
the run is not in the file.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `acc` | `rod.meta` | derived | Accession of the SRA run: the first part of the contig name (`DRR000601` in `DRR000601_6173_circle_6173`). |
| `organism` | `rod.meta` | file | Organism, or metagenome, the library was made from. |
| `assay_type` | `rod.meta`, `rod.df` | file | Library strategy, e.g. `WGS` (genomic DNA) or `RNA-Seq` (transcribed RNA). |
| `librarysource` | `rod.meta` | file | Library source, e.g. `GENOMIC`, `METAGENOMIC`. |
| `biosample` | `rod.meta` | file | BioSample accession of the run. |
| `bioproject` | `rod.meta` | file | BioProject accession of the run. |

## RNA-Seq enrichment scores

Computed by `02A_feature_tests.Rmd` and `03_enriched_contigs.Rmd`.

- **Groups.** Contigs of the two most common assay types (WGS, RNA-Seq) are
  split into four groups by assay type and `circular`.
- **Test.** Each of 13 features is binned. The WGS circular contigs are the
  reference: they are down-sampled 1,000 times to the size of a group, which
  gives the expected count (mean) and its standard deviation for every bin.
- **Bin scores.** Z-score = (observed - expected) / standard deviation;
  fold-change = observed / expected.
- **Contig scores.** A contig takes the score of its bin in each feature, and
  the scores are summed over the 14 features.
- **Features tested:** `length`, `mfe`, `mfent`, `mfeZ`, `gc_content`,
  `orf_length`, `orf_pct`, `bKeto`, `uKeto`, `entropy`, `meanEntropy20`,
  `pRod`, `nBranch`, `length_mod3`.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `group` | `rod.score` | derived | Group the contig was tested in: assay type followed by `circular` or `other`, e.g. `RNA-Seq circular`. `NA` for the other assay types, which have no scores. |
| `z` | `rod.score` | derived | Sum of the bin Z-scores. |
| `z_norm` | `rod.score` | derived | Sum of the bin Z-scores, each limited to 0-10. |
| `fc` | `rod.score` | derived | Sum of the bin fold-changes. 14 for a WGS circular contig (observed = expected in every bin). |
| `fc_norm` | `rod.score` | derived | Sum of the bin fold-changes, each limited to 0-10. |
| `sig` | `rod.score` | derived | Sum of the bin fold-changes limited to 0-10, counting only bins with a Z-score above 2. |

## UMAP

Computed by `06_umap.Rmd`: a 2D map (UMAP) of all contigs by 33 features:
`length`, `length_mod3`, `gc_content`, `entropy`, `minEntropy20`,
`meanEntropy20`, `circular`, `mfeZ`, `orf_pct`, `n_orfs`, `orf_strand`, and
eleven of the fold of each strand: `mfent`, `pS`, `pH`, `pI`, `pB`, `pF`,
`pRod`, `nBranch`, `pLongestRod`, `bKeto`, `uKeto` and the same with `_rc`
for the reverse complement. Each is transformed and scaled to mean 0 and
standard deviation 1. The transformations are `log10` for `length` and
`asinh` (a signed logarithm) for `mfeZ`; every other feature is used as it
is. The map and its clusters are made from the ROD contigs alone (seed
31337, Leiden resolution 3.5, which gives 49 clusters in snapshot `v261126`,
8 of them RNA-Seq-enriched; the resolution of the earlier maps, 2.5, gives
38). Every ROD contig has a position: a contig without a base pair on a
strand has 0 for the structure features of that strand (the 9,931 without
one on the reverse strand are a cluster of their own, 48). Earlier snapshots
had other maps, which this one does not follow: `v261120` (14 features of the forward
strand, 27 clusters), `v261122` (20 features, 50 clusters), `v261123` (21
features, with `mfe` and `mfe_rc` in place of `mfent`, `mfent_rc` and
`meanEntropy20`, and the old ORF calls; 47 clusters) and `v261124` and
`v261124b` (22 features: `keto_skew` of the whole contig, and of each
strand `mfent`, `pct_backbone`, `n_branch_points`, `pct_bound` and
`pct_loop`; resolution 2.5, 45 clusters, 7 enriched; a contig without a
base pair on either strand had no position). The map is fitted on a
random sample of 500,000 contigs; the other contigs are placed on it by their
nearest neighbours in the sample. Contigs with similar features lie close
together; the axes have no meaning of their own.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `umap_x`, `umap_y` | `rod.umap` | derived | Position of the contig on the map. |
| `umap_cluster` | `rod.umap` | derived | Cluster of the contig, 1 the largest: Leiden clusters of the neighbour graph of the sample, passed on to the other contigs by the most common cluster of their nearest neighbours. |

The control sequences take no part in the map. Each is placed on it
afterwards, by its nearest neighbours in the sample, when it has all 33
features; its position and cluster are in
`control.umap`, with the columns of `rod.umap` (`NA` for a control without
`mfeZ`, which lacks a feature). The fitted map is saved as
`umap_model.uwot`, with the scaling of the features in `umap_features`.
`umap3d` holds a second UMAP, in three dimensions, of a sample of the contigs
and controls (`source`, `source_row`, `x`, `y`, `z`).

The table `umap_clusters` has one row per cluster: its number of contigs
(`contigs`), its centre on the map (`x`, `y`), the share of circular contigs
(`pct_circular`), its RNA-Seq and WGS contigs (`n_test`, `n_reference`), the
RNA-Seq contigs expected from the WGS contigs (`expected`), the fold-change
observed / expected and its Z-score (`fold`, `z`), the same over circular
contigs only (`fold_circular`, `z_circular`), and `enriched` (Z above 2 and
fold-change of 2 or more). For each control category it has the number of
controls in the cluster and their fold-change and Z-score against the WGS
circular contigs (`n_vle`, `fold_vle`, `z_vle`, `n_decoy`, `fold_decoy`,
`z_decoy`).

## RNA-cluster UMAP

Computed by `07_umap_rna_cluster.Rmd`: a second map (UMAP2) of the contigs
in the RNA-Seq-enriched clusters of the first map alone (clusters with a
Z-score above 2 and a fold-change of 2 or more), from the same 33 features,
scaled again within these contigs. A contig that is not in one of those
clusters has `NA` in every column. The controls of the same clusters are
placed on the map afterwards (`control.umap_rna`); they take no part in it. In
snapshot `v261126` the map holds 2,162,644 contigs of 8 clusters and has
36 clusters of its own, 17 of them RNA-Seq-enriched (in `v261124` and
`v261124b`: 1,673,544 contigs of 7 clusters, 41 clusters).

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `rna_umap_x`, `rna_umap_y` | `rod.umap_rna` | derived | Position of the contig on the RNA-cluster UMAP. |
| `rna_umap_cluster` | `rod.umap_rna` | derived | Cluster of the contig on that map, 1 the largest (Leiden clusters of its neighbour graph). |
| `RNA_cluster` | `rod.umap_rna` | derived | Label of that cluster: its number and, where it has one, its enriched element, the kind of control with the highest fold-change in the cluster against the contigs of the map (10 controls or more, fold of 2 or more, Z above 2), e.g. `22 Ribozyvirus`; the number alone for a cluster without one. |

The table `umap_rna_clusters` has one row per cluster of this map: its
contigs, its centre (`x`, `y`), the RNA-Seq test against all ROD contigs
(`n_test`, `n_reference`, `fold`, `z`, `enriched`) and within the map
(`fold_within`), the SRA runs of its RNA-Seq contigs and the share of the
largest run and BioProject (`n_runs`, `top_run_pct`, `top_project_pct`), the
cluster of the first map most of its contigs come from (`origin`,
`origin_pct`), its controls (`n_vle`, `n_decoy`), its enriched element
(`element`, `element_n`, `element_fold`, `elements`) and its label
(`RNA_cluster`).

## ORF clusters

Made by `08_orf_clusters.Rmd` from the ORF clusters of `data/cluster/`: every
ORF of 60 codons or more of the ROD contigs and of the controls, clustered
by protein sequence with DIAMOND at 70% identity (`cluster70`), and the
centroids of those clusters at 30% identity (`cluster30`), with a member and
its centroid covering 80% of each other. A 70% cluster belongs to one 30%
cluster. Until snapshot `v261124` the tables had `select_cluster` in the
place of these: a summary of clusters made outside the project from an
earlier set of ORF calls.

**One row per ORF** (`orf.df`), the ROD ORFs first, in the order of the ORF
files; the ORFs of a contig are next to each other.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `category` | `orf.df` | file | `ROD`, `DECOY` or `VLE`. |
| `contig_row` | `orf.df` | derived | Row of the contig of the ORF in the `rod` tables (`ROD`) or in the `control` tables. |
| `orf_nr` | `orf.df` | file | Number of the ORF on its contig, as in the ORF file (`orf3`). |
| `orf_start`, `orf_end`, `strand`, `frame`, `orf_code` | `orf.df` | file | Position of the ORF and its ORF code, as in the ORF file. |
| `aa_length` | `orf.df` | derived | Length of the protein. |
| `cluster70`, `cluster30` | `orf.df` | derived | Its clusters: rows of `orf.cluster70` and of `orf.cluster30`. |
| `is_centroid` | `orf.df` | derived | Is the ORF the centroid of its 70% cluster? |
| `identity`, `member_cover`, `centroid_cover` | `orf.df` | derived | Best local alignment of the ORF to the centroid of its 70% cluster: identical positions (%), and the part of the ORF and of the centroid in the alignment (%). |

**One row per cluster** (`orf.cluster70`, `orf.cluster30`); the row number is
the number of the cluster, 1 the cluster with the most ORFs. A contig of a
cluster is a ROD contig with at least one ORF of the cluster.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `centroid` | both | derived | Row of the centroid ORF in `orf.df`. |
| `centroid_category`, `centroid_contig_row`, `centroid_orf_nr`, `centroid_aa_length` | both | derived | The centroid ORF, as in `orf.df`: the set of its contig (`ROD`, `DECOY`, `VLE`), the row of the contig in the tables of that set, the number of the ORF on the contig and the length of its protein. The name of the centroid (`DRR000346_804025_orf1`) is put together from them where it is shown. |
| `primary_orf` | both | derived | Is the centroid ORF the longest ORF of its contig (the ORF of `orf_length` in `rod.orf`; the first of them when two have that length)? The explorer lists the clusters of a primary ORF by default: a contig mostly has other, shorter ORFs, each the centroid of a cluster of its own. |
| `cluster30` | `orf.cluster70` | derived | The 30% cluster the 70% cluster belongs to. |
| `clusters70` | `orf.cluster30` | derived | Number of 70% clusters in the 30% cluster. |
| `members` | both | derived | ORFs in the cluster. |
| `n_rod`, `n_decoy`, `n_vle` | both | derived | Of them: ORFs of ROD contigs, of decoy controls and of VLE controls. |
| `contigs` | both | derived | ROD contigs of the cluster. |
| `n_test`, `n_reference`, `expected` | both | derived | Its RNA-Seq and WGS contigs, and the RNA-Seq contigs expected from its WGS contigs (a down-sample of the WGS contigs with an ORF to the number of RNA-Seq contigs with one). |
| `fold`, `z` | both | derived | Fold-change observed / expected and its Z-score, as for the UMAP clusters. `fold` is infinite for a cluster with RNA-Seq and without WGS contigs; both are missing where they cannot be worked out. |
| `fold_circular`, `z_circular` | both | derived | The same over circular contigs only. |
| `enriched` | both | derived | Z above 2 (or no WGS contig), a fold-change of 2 or more, and 10 or more RNA-Seq and WGS contigs. |
| `meanLength`, `meanGC`, `meanMFE` | both | derived | Mean of `length`, `gc_content` and `mfe` over its contigs. |
| `sdLength` | both | derived | Standard deviation of `length` over its contigs; missing for a cluster of one contig. |
| `meanRod`, `sdRod` | both | derived | Mean and standard deviation of `pRod` (forward strand) over its contigs. |
| `n_libraries`, `n_bioprojects` | both | derived | Distinct SRA runs and BioProjects of its contigs. |
| `percent_RNA` | both | derived | Percentage of its libraries whose assay type is `RNA-Seq`. |
| `organism`, `organism_pct` | both | derived | The most common organism among its contigs, and the percentage of the contigs from it (until 2026-11-24 in the evening: among its libraries). |
| `organisms` | both | derived | The organisms of its contigs with the percentage of the contigs from each, the most common first: `Homo sapiens(97%); gut metagenome(3%)`. The first five are named, the rest counted (`12 others(8%)`); a contig whose organism is missing would be `unknown` (in ROD v1.1 none is: 16,253 contigs have the text `NA` as their organism in the metadata, and are listed as `NA`); a percentage under 0.5 is `<1`. |
| `umap_cluster`, `umap_cluster_pct` | both | derived | The cluster of the UMAP with the most of its contigs, and the percentage of its contigs there. |
| `umap_enriched`, `umap_enriched_pct`, `umap_enriched_fold` | both | derived | The cluster of the UMAP in which its contigs are most over-represented: the one least likely to hold so many of them by chance, if they were drawn from all contigs with an ORF and a position on the map (hypergeometric). With it, the percentage of its contigs there, and that share over the share of the map cluster among all such contigs. |
| `rna_umap_cluster`, `rna_umap_cluster_pct` | both | derived | The same for the RNA-cluster UMAP (the percentage is of all its contigs). |
| `model`, `model_evalue`, `model_bitscore` | both | derived | The protein model of the cluster, with the E-value and bit score of the match: for a 70% cluster the model that its centroid matches significantly (`hmmsearch`, E-value of 0.01 or less; `AmRdRp`, `dAg`, `Oblin1`, `Oblin2`), for a 30% cluster the model that most of its matching 70% clusters have, with their best match. |
| `model_clusters70` | `orf.cluster30` | derived | Number of its 70% clusters with a model. |

**One row per contig** (`rod.orfcluster`, `control.orfcluster`): the ORF
annotation of the contig, which the explorer shows.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `orf_first` | `rod.orfcluster` | derived | Row in `orf.df` of the first ORF of the contig: its ORFs are that row and the `n_orfs - 1` after it. Missing for a contig without an ORF of 60 codons. |
| `orf_cluster70`, `orf_cluster30` | `rod.orfcluster` | derived | The clusters of the longest ORF of the contig. |
| `orf_cluster_members` | `rod.orfcluster` | derived | ORFs in that 30% cluster. |
| `orf_cluster_fold`, `orf_cluster_enriched` | `rod.orfcluster` | derived | RNA-Seq fold-change of that 30% cluster, and whether it is enriched. |
| `orf_model` | `rod.orfcluster` | derived | The protein models of the ORFs of the contig (the model of the 70% cluster of each ORF). A contig with ORFs of two models has both, with a `+` between them (`Oblin1+Oblin2`). |

