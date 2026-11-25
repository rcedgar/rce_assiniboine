# Circular RNA Longest-Stem Parser & Residue Annotator

## Overview

This algorithm parses a 3-symbol dot-bracket string (`(`, `)`, `.`) representing a **circular RNA secondary structure**, identifies the longest composite stem, and produces residue-level structural annotations for both paired and unpaired nucleotides.

Unlike naive linear parsers, this algorithm accounts for:

1. **Circular topology:** The 5' and 3' ends are covalently linked, allowing helical stems to wrap seamlessly across the `N - 1` $\rightarrow$ `0` index boundary.
2. **Composite helices (bulges & internal loops):** Consecutive base-paired helices interrupted by small unpaired regions (up to `max_unpaired_gap` nucleotides per strand) are merged into a single continuous stem, while multi-branch junctions act as strict terminators.
3. **Full residue classification:** Every nucleotide—paired or unpaired—is assigned a structural role and flagged if it belongs to the longest stem.

---

## Algorithm Workflow

### 1. Stack-Based Base-Pair Extraction

The input string of length $N$ is scanned linearly from index `0` to `N - 1` using a LIFO stack:

* `(` pushes its index onto the stack.
* `)` pops the most recent `(` index `j` and records a bidirectional mapping (`pair_map[i] = j` and `pair_map[j] = i`).
* `.` is skipped during initial pairing.
* Unmatched brackets or invalid characters raise a `ValueError`.

### 2. Circular Modulo Indexing & Bidirectional Stem Traversal

Because a circular sequence can be linearized at an arbitrary cut site—potentially splitting a single stem across indices `N - 1` and `0`—all index increments and decrements use modulo arithmetic (`(idx ± step) % N`).

For each unvisited canonical pair $(i, j)$:

1. **Outward Walk (`prev_outward_pair`):** Steps outward toward the exterior/junction (`i - step_5`, `j + step_3`) until it hits either a gap exceeding `max_unpaired_gap` or a branch point, locating the true outermost closing pair of the stem.
2. **Inward Walk (`next_inward_step`):** Steps inward toward the apical loop (`i + step_5`, `j - step_3`), collecting all stacked base pairs and recording any intermediate unpaired indices bridged along both strands.
3. **Symmetric Deduplication:** Both `(i, j)` and `(j, i)` are added to a `visited` set so stems crossing the origin boundary are never traversed twice.

### 3. Branch-Point Detection vs. Internal Loops

When stepping from base pair $(i, j)$ to the next paired nucleotides encountered (`ni` on the 5' strand and `nj` on the 3' strand within `max_unpaired_gap`):

* **If `pair_map[ni] == nj`:** The intervening dots form a bulge or internal loop. The helices are merged, and the intervening dots are recorded as `internal_unpaired` residues of that stem.
* **If `pair_map[ni] != nj`:** The stem has reached a **multi-branch junction** (the two strands diverge into separate hairpins). Traversal stops immediately so distinct branches are not merged.

### 4. Ranking & Residue Role Classification

All discovered stems are ranked by their total count of paired bases (`len(stem["pairs"])`). Each nucleotide `0 .. N - 1` is then classified into one of four structural roles:

| Role | Condition | Included in `in_longest_stem`? |
| --- | --- | --- |
| **`paired_stem`** | Nucleotide has a base-pairing partner. | `True` if part of the winning stem's base pairs. |
| **`unpaired_bulge_or_internal_loop`** | Unpaired dot bridged between two helices of a composite stem. | `True` if bridged within the winning stem. |
| **`unpaired_hairpin_loop`** | Unpaired dots enclosed inside the innermost pair of a stem. | Controlled by `include_hairpin_loop` (default `False`). |
| **`unpaired_junction_or_exterior`** | Unpaired dots in multi-branch junctions or exterior regions. | Always `False`. |

---

## Parameters & Return Schema

### Input Parameters

* **`dot_bracket`** (`str`): Secondary structure string composed of `(`, `)`, and `.`.
* **`max_unpaired_gap`** (`int`, default `2`): Maximum number of consecutive unpaired nucleotides allowed on either strand when bridging adjacent helices.
* `0`: Strict contiguous helices only (any bulge or internal loop splits the stem).
* $\ge 1$: Merges helices separated by bulges or symmetric/asymmetric internal loops up to that size.


* **`include_hairpin_loop`** (`bool`, default `False`): Whether apical hairpin loop dots capping the innermost pair of the longest stem should be marked `in_longest_stem = True`.

### Return Dictionary

```json
{
  "longest_stem_length": 5,
  "longest_stem_pairs": [[0, 16], [1, 15], [4, 12], [5, 11], [6, 10]],
  "longest_stem_unpaired": [2, 3, 13, 14],
  "annotations": [
    {
      "index": 0,
      "symbol": "(",
      "is_paired": true,
      "partner": 16,
      "role": "paired_stem",
      "in_longest_stem": true
    }
  ]
}

```

---

## Complexity & Edge Cases

* **Time Complexity:** $O(N \cdot G^2)$ where $N$ is sequence length and $G$ is `max_unpaired_gap` (effectively $O(N)$ for small constant $G$).
* **Space Complexity:** $O(N)$ to store `pair_map`, `visited`, and the output annotation array.
* **Fully Unpaired Structures:** Strings containing only `.` safely return `longest_stem_length = 0` and classify all residues as `unpaired_exterior`.
* **Tie-Breaking:** If multiple stems share the maximum base-pair count, `max()` selects the first stem encountered in canonical 5' order.