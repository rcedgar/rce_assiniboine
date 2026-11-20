# Features of ROD v1.1

Every feature of the 11,308,939 contigs of ROD v1.1 used in this project:
the values read from the original data files and the values calculated from
them by the `rodFold` package and the notebooks of this repository.

- **Origin** says which is which: `file` for a value read unchanged from a
  data file, `derived` for a value calculated in this project.
- **One row per contig.** Every table has the same rows in the same order
  (the order of `ROD_v1.1.fa`), so row `i` is the same contig in each table.
  The contig names are in `rod.df$contig`.
- **Table** is the `.RData` file that holds the feature (`rod.df` is saved as
  `rod.df.RData`, and so on). Several features are also copied into `rod.df`.
- **Percentages** run from 0 to 100 and are of the contig length unless
  stated otherwise.

| Data file | Holds |
|:----------|:------|
| `ROD_v1.1.fa` | Contig sequences, with the k-mer coverage and assembly-graph edges in the header. |
| `ROD_v1.1_dotbracket.tsv` | RNAfold prediction of each contig: minimum free energy and structure. |
| `ROD_v1.1.reorf_n60.fa` | Open reading frames called on the contigs, one record per ORF. |
| `ROD_v1.1_SRARunInfo.csv` | Metadata of the SRA run each contig was assembled from. |
| `rce_cluster_summary.fev` | ORF clusters of the contigs (ORFs clustered at 70% amino acid identity): one row per cluster with its number of members. |
| `rce_cluster_2026-11-16/` | One folder per ORF cluster with 10 or more members: its member contigs (`circles.fa`) and their nucleotide circular clusters (`circuclust.tsv`). |

## Contig and sequence

Read from `ROD_v1.1.fa`. Each record is a BCALM2 unitig with a header such as
`>DRR000601_6173_circle_6173 ka:f:9.342 L:-:6173:- L:+:6173:+`.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `contig` | `rod.df` | file | Contig name, `<SRA run>_<unitig id>_circle_<unitig id>`: the header up to the first space. |
| `header` | `rod.fa` | file | Full header line without the `>`: contig name, k-mer coverage and assembly-graph edges. |
| `sequence` | `rod.fa` | file | Nucleotide sequence of the contig, as written in the FASTA. |
| `length` | `rod.df` | file | Contig length: number of bases of the sequence (nt). |
| `kmer_coverage` | `rod.df` | file | k-mer coverage (`ka`): the number of the header field `ka:f:`, so `ka:f:9.342` is 9.342. BCALM2 reports it for every unitig: the mean abundance of its k-mers in the reads of the library, i.e. how deeply the contig was sequenced. |

## Sequence composition

Computed from the nucleotide sequence by `read_rod_fasta()`. `U` counts as
`T`, case is ignored, and other letters (e.g. `N`) are left out.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `gc_content` | `rod.df` | derived | G + C bases as a percentage of the contig length. |
| `keto_skew` | `rod.df` | derived | Absolute keto skew, `abs(GT - AC) / (GT + AC)` with `GT` and `AC` the numbers of keto (G, T) and amino (A, C) bases. 0 when both strands hold the same number of keto bases, 1 when one strand holds them all. The same for a contig and its reverse complement. |
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

The RNAfold prediction of each contig, read from `ROD_v1.1_dotbracket.tsv`,
and the stability measures calculated from it.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `mfe` | `rod.fold`, `rod.df` | file | Minimum free energy (MFE) of the predicted fold (kcal/mol). More negative is more stable. |
| `dot_bracket` | `rod.dotbracket` | file | Predicted secondary structure in dot-bracket notation, one character per base: `(` and `)` for the two bases of a pair, `.` for an unpaired base. |
| `mfent` | `rod.fold` | derived | Minimum free energy per nucleotide, `mfe / length` (kcal/mol/nt). MFE grows with length; this is the length-independent measure of fold stability. |
| `mfentgc` | not saved | derived | `log(-mfent / gc_content)` (natural log): fold stability relative to GC content. Calculated for the figures of `02B_explore_features.Rmd` only. |

## Structure composition

Computed from the dot-bracket structure by `parse_dot_bracket()`. Every base
is one of three kinds, so `n_bound + n_loop + n_unbound = length`. A hairpin
is one top-level opening and closing of brackets; branches nested inside it
belong to it.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `n_bound`, `pct_bound` | `rod.fold` | derived | Paired bases, `(` or `)`: number, and percentage of the contig length. |
| `n_loop`, `pct_loop` | `rod.fold` | derived | Unpaired bases enclosed by at least one base pair (inside a hairpin): number, and percentage. |
| `n_unbound`, `pct_unbound` | `rod.fold` | derived | Unpaired bases outside every base pair (before, after or between hairpins): number, and percentage. |
| `n_hairpins` | `rod.fold` | derived | Number of hairpins. |
| `longest_hairpin`, `pct_longest_hairpin` | `rod.fold` | derived | Length of the longest hairpin, from its opening to its closing bracket (nt), and as a percentage of the contig length. |

## Rod-likeness

Computed from the dot-bracket structure by `dot_bracket_branching()`. The
structure is read as a circular molecule, as a tree of loops joined by base
pairs.

- A loop with three or more base pairs at its edge is a **branch point**.
- The bases outside every pair count as one more loop, so two hairpins side by
  side are one stem that runs across position 1.
- An **unbranched stem** is a run of base pairs joined only by bulges and
  interior loops, from one end loop or branch point to the next.
- The **backbone** is the path from one end loop to another that holds the
  most base pairs; what hangs off it is a **side branch**.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `pct_backbone` | `rod.fold` | derived | Backbone fraction: base pairs on the backbone as a percentage of all base pairs. 100 for a rod. |
| `n_branch_points` | `rod.fold` | derived | Number of branch points. 0 for a rod. |
| `pct_branch` | `rod.fold` | derived | Bases in side branches as a percentage of the contig length. |
| `unbranched_length`, `pct_unbranched` | `rod.fold` | derived | Bases in the longest unbranched stem (its paired bases, the loops between them and its end loops): number, and percentage of the contig length. |
| `unbranched_bound`, `pct_unbranched_bound` | `rod.fold` | derived | Paired bases in that stem: number, and percentage of the contig length. |

## Open reading frames

Read from `ROD_v1.1.reorf_n60.fa`, which has one record per ORF of 180 nt (60
codons) or longer, with a header such as
`>ROD|DRR235453_664787_circle_664787 orf1 code1 3-329/327(+3)Cem`: contig,
ORF number, genetic code table, `start-end`, contig length, `(strand, frame)`
and an optional ORF code. The protein sequences are not kept.

The tables hold the **longest ORF** of each contig: the record with the
largest length, and on ties the first in the file. A contig without a record
has no ORF of 180 nt or longer.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `orf_length` | `rod.orf`, `rod.df` | file | Length of the longest ORF, `end - start + 1` (nt); 0 if the contig has no ORF. |
| `orf_strand` | `rod.orf`, `rod.df` | file | Strand of the longest ORF, `+` or `-`; `NA` if none. |
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
  the scores are summed over the 13 features.
- **Features tested:** `length`, `mfe`, `mfent`, `gc_content`, `orf_length`,
  `orf_pct`, `keto_skew`, `entropy`, `meanEntropy20`, `deltaEntropy20`,
  `pct_longest_hairpin`, `pct_backbone`, `n_branch_points`.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `group` | `rod.score` | derived | Group the contig was tested in: assay type followed by `circular` or `other`, e.g. `RNA-Seq circular`. `NA` for the other assay types, which have no scores. |
| `z` | `rod.score` | derived | Sum of the bin Z-scores. |
| `z_norm` | `rod.score` | derived | Sum of the bin Z-scores, each limited to 0-10. |
| `fc` | `rod.score` | derived | Sum of the bin fold-changes. 13 for a WGS circular contig (observed = expected in every bin). |
| `fc_norm` | `rod.score` | derived | Sum of the bin fold-changes, each limited to 0-10. |
| `sig` | `rod.score` | derived | Sum of the bin fold-changes limited to 0-10, counting only bins with a Z-score above 2. |

## ORF cluster summaries

One row per ORF cluster, not per contig. The clusters (ORFs clustered at 70%
amino acid identity) are read from `rce_cluster_summary.fev` and the folder
`rce_cluster_2026-11-16/`. `04_cluster_analysis.Rmd` keeps the clusters with
10 or more members whose centroid is a circular contig of an RNA-Seq library,
and summarises their members that are ROD contigs.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `cluster_nr` | `select_cluster` | file | Number of the ORF cluster. |
| `members` | `select_cluster` | file | Number of members of the ORF cluster. |
| `groups` | `select_cluster` | file | Sequence sets of the members and the number from each, e.g. `ROD:1062;`. |
| `centroid` | `select_cluster` | file | Contig of the cluster centroid: the centroid of its first nucleotide cluster. |
| `nt_clusters` | `select_cluster` | derived | Number of nucleotide circular clusters in the ORF cluster. |
| `n_rod` | `select_cluster` | derived | Members that are ROD contigs. |
| `meanLength`, `meanGC`, `meanMFE`, `meanHairpins` | `select_cluster` | derived | Mean of `length`, `gc_content`, `mfe` and `n_hairpins` over the members. |
| `n_libraries`, `n_bioprojects` | `select_cluster` | derived | Number of distinct SRA runs and BioProjects of the members. |
| `percent_RNA` | `select_cluster` | derived | Percentage of the libraries of the members whose assay type is `RNA-Seq`. |
| `centroid_organism` | `select_cluster` | derived | Organism of the library of the centroid contig. |
| `organisms` | `select_cluster` | derived | Organisms of the libraries of the members with the percentage of libraries from each, most common first. |
