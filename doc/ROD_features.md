# Features of ROD v1.1

## Sequence composition

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `gc_content` | `rod.df` | derived | G + C bases as a percentage of the contig length. |
| `keto_skew` | `rod.df` | derived | Absolute keto skew, `abs(GT - AC) / (GT + AC)` with `GT` and `AC` the numbers of keto (G, T) and amino (A, C) bases. 0 when both strands hold the same number of keto bases, 1 when one strand holds them all. The same for a contig and its reverse complement. |
| `entropy` | `rod.entropy`, `rod.df` | derived | Shannon entropy of the base composition of the whole contig, `-sum(p * log2(p))` over A, C, G, T (bits). 2 when the four bases are equally common, 0 for a homopolymer. |
| `minEntropy20` | `rod.entropy` | derived | Lowest Shannon entropy of any 20 nt window of the contig (bits): the least complex 20 nt stretch. Windows start at every base and do not wrap around the end of the contig. |
| `meanEntropy20` | `rod.entropy` | derived | Mean Shannon entropy of all 20 nt windows (bits). |
| `deltaEntropy20` | `rod.entropy` | derived | `minEntropy20 - entropy` (bits). 0 or below; strongly negative when a contig that uses all four bases has a low-complexity stretch. |

## Structure composition

Computed from dot-bracket string.
Every base is one of three kinds, so `n_bound + n_loop + n_unbound = length`.
A hairpin is one top-level opening and closing of brackets; branches nested inside it belong to it.

| Feature | Table | Origin | Description |
|:--------|:------|:-------|:------------|
| `n_bound`, `pct_bound` | `rod.fold` | derived | Paired bases, `(` or `)`: number, and percentage of the contig length. |
| `n_loop`, `pct_loop` | `rod.fold` | derived | Unpaired bases enclosed by at least one base pair (inside a hairpin): number, and percentage. |
| `n_unbound`, `pct_unbound` | `rod.fold` | derived | Unpaired bases outside every base pair (before, after or between hairpins): number, and percentage. |
| `n_hairpins` | `rod.fold` | derived | Number of hairpins. |
| `longest_hairpin`, `pct_longest_hairpin` | `rod.fold` | derived | Length of the longest hairpin, from its opening to its closing bracket (nt), and as a percentage of the contig length. |

## Rod-likeness

Computed from the dot-bracket string. The structure is read as a circular molecule,
as a tree of loops joined by base pairs.

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
