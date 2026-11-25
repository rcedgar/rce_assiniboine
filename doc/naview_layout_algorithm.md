# NAView layout

Layout of an unknotted RNA secondary structure as a modified radial drawing.
The method is the one in Bruccoleri and Heinrich, *Computer Applications in
the Biosciences* 4:167–173 (1988), as implemented by `naview/naview.c`.

The programming language is not fixed. This document is the contract: the
records below, and the geometry that fills the coordinates.

## Input

One string of length `n ≥ 1`, and a flag `circular` (default false). Every
character of the string is `(`, `.`, or `)`.

| Character | Meaning |
|:----------|:--------|
| `(` | 5' base of a pair |
| `)` | 3' base of that pair |
| `.` | unpaired base |

Brackets are balanced and properly nested. The alphabet cannot express a
pseudoknot. A string that uses any other character, or whose brackets do not
match, is rejected.

Bases are numbered `0 .. n-1` in 5' to 3' order. `start` and `end` on every
record are inclusive. A 1-based view of the same base is index `+ 1`. When
`circular` is false, base `n-1` and base `0` have no backbone bond. When
`circular` is true, that bond is real. The brackets stay a balanced linear
string either way: the flag does not add pairs, and a helix still cannot
cross the numbering origin. A caller who wants the origin inside a loop
rotates the string until the brackets balance, then sets the flag.

## Output

Coordinates, loops, and segments are always returned. Loops and segments are
determined by the dot-bracket and by `circular`. Which loop is drawn at the
origin affects coordinates only.

Two further strings are optional:

| Option | String | What it is |
|:-------|:-------|:-----------|
| NAView graph | `naview_graph` | The loop tree, labeled and in backbone order. |
| RAG | `rag_tree` and `rag_dual` | The same tree after the RNA-As-Graphs contractions, then the dual of that graph. |

The grammars are in [Graph notation](#graph-notation). Neither string depends
on the drawing root or on the nucleotide coordinate frame. Each graph may
also return a picture of itself, in [Graph layout](#graph-layout): one point
per vertex, and one line record `(x0, y0, x1, y1)` per drawn stroke.

### Coordinates

For each base `i`, a pair `(x, y)`.

Scale, origin, rotation, and which way is up are arbitrary. A correct layout
may be translated, rotated, uniformly scaled, or reflected, and it is still
the same result. Absolute numbers cannot be compared across implementations.
What is fixed is the drawing: stems are straight, a base pair has the same
length as a target backbone step, loops are circular except where a strand is
extruded, and stem directions come from the circle diagram below.

A reference implementation can use the historical picture frame if it wants
one: backbone length near 1, root loop centered at the origin, `+x` to the
right, `+y` up, angles increasing counterclockwise. Callers must treat that
frame as a choice.

### Loops

One record per loop.

| Field | Meaning |
|:------|:--------|
| `loop_id` | `0` for the first loop created, then `1, 2, …` in the preorder of [Loop tree](#loop-tree). Loop `0` is the exterior when `circular` is false, and the closing loop when `circular` is true. |
| `kind` | `exterior`, `circle`, `hairpin`, `bulge`, `interior`, or `multibranch` |
| `start`, `end` | Inclusive span. See below. |
| `parent_loop_id` | The loop on the other side of the closing helix. Absent for the exterior. |
| `closing_helix_id` | Helix that closes this loop. Absent for the exterior. |
| `child_helix_ids` | Helices that leave this loop toward nested loops, in backbone order. |

Span:

- Exterior (`circular` false): `start = 0`, `end = n-1`. This span is the
  whole chain because the exterior lies outside every pair. The bases that
  belong to the exterior are only its unpaired segments, and those segments
  may be disjoint. The 5' tail and the 3' tail are separate segments; the
  layout joins them by a phantom slot that is not a base.
- Closing loop (`circular` true): no closing pair and no parent. `start` and
  `end` name the unpaired segment that contains the bond from `n-1` to `0`.
  When that bond has unpaired bases on both sides, `start > end`: the bases
  are `start .. n-1` and `0 .. end`. When the bond has no unpaired base,
  `start` and `end` are absent. Nested helices are listed only in
  `child_helix_ids`.
- Every other loop: `start` and `end` are the two bases of the closing pair
  (the inner pair of `closing_helix_id`). They pair with each other. The open
  range `(start, end)` contains that loop's unpaired segments and everything
  nested inside it.

Kind, from the helices that touch the loop:

| Kind | When |
|:-----|:-----|
| `exterior` | `circular` is false, and this is the unique loop outside every pair. The kind stays `exterior` however many stems it has, including zero. |
| `circle` | `circular` is true and the string has no pair. One loop, spanning every base. |
| `hairpin` | One helix. On the closing loop this is the outer pair of that helix; the unpaired bases lie outside that pair and contain the bond from `n-1` to `0`. |
| `bulge` | Two helices, and exactly one of the two gaps between them contains no base. |
| `interior` | Two helices, and the two gaps are not a one-empty bulge. Both gaps may contain bases. Both gaps empty occurs on the closing loop when two top-level helices meet on both sides with no unpaired base; the helix finder leaves them as two helices because they are not a nested stack. |
| `multibranch` | Three or more helices. |

A stack of consecutive pairs is one helix, so a loop is never emitted between
stacked pairs. A hairpin, bulge, or interior side with no unpaired base is
omitted from the segment list; the loop record is still emitted. `((()))` is
one helix, an exterior, and a hairpin whose unpaired segment is absent.

### Segments

One record per contiguous run a post-processor might color. Each segment has
a single inclusive interval.

| Field | Meaning |
|:------|:--------|
| `segment_id` | `0, 1, …` in increasing order along the backbone, starting at base `0` |
| `kind` | `helix-5p`, `helix-3p`, `hairpin`, `bulge`, `interior`, `multibranch`, `exterior`, or `circle` |
| `start`, `end` | Inclusive. `start ≤ end` is one linear run. `start > end` is a run that contains the bond from `n-1` to `0`: the bases are `start .. n-1` and `0 .. end`. That form is used only when `circular` is true. |
| `helix_id` | Set on both strands of a helix. The same id ties the two strands together. |
| `loop_id` | Set on every unpaired segment. Several segments can share one loop. |
| `side` | `5p` or `3p` for a helix strand and for each side of a bulge or interior loop. Absent for hairpin, multibranch, and exterior. |

| Kind | Interval |
|:-----|:---------|
| `helix-5p` | 5' strand of one helix (the lower indices). |
| `helix-3p` | 3' strand of that helix (the higher indices). |
| `hairpin` | The unpaired bases of a hairpin, strictly inside its closing pair. |
| `bulge` | The one non-empty unpaired side of a bulge. |
| `interior` | One unpaired side of an interior loop. The loop has two of these, distinguished by `side`. |
| `multibranch` | One unpaired gap between consecutive helices of a multibranch loop. |
| `exterior` | One unpaired piece of the exterior when `circular` is false: a 5' tail, a 3' tail, or a stretch between two top-level helices. |
| `circle` | The single unpaired run when `circular` is true and the string has no pair. |

Empty gaps are not segments. When `circular` is true the 5' tail and the 3' tail of the closing loop are one segment, because the bond between them is real; `start > end` when both tails contain a base. The two strands of a helix are two segments.
Both sides of an interior loop are two segments. A multibranch loop contributes
one segment per gap that contains a base.

A helix has no record of its own. Its identity sits on the two strand
segments: shared `helix_id`, the two intervals, and the enclosing span
`[start of helix-5p, end of helix-3p]`. The inner loop of that helix is the
loop whose `closing_helix_id` equals the helix. The outer loop is that inner
loop's `parent_loop_id`. Number helices `0, 1, …` by increasing 5' start.

### Nesting

Two different interval families exist.

**Segments partition the sequence.** Every base belongs to exactly one
segment. Segment intervals meet end to end, do not overlap, and do not
contain one another. Coloring bases from the segment list needs no paint
order. Paired bases take a helix-strand color. Unpaired bases take the color
of their loop segment. To color a whole stem, color both segments that share
its `helix_id`. To color a whole loop's own bases, color the unpaired
segments with that `loop_id` together with the two bases of its closing pair.

**Loop spans and helix spans nest.** They form a tree. Where two of these
spans overlap, one contains the other. The exterior span is the whole chain
and contains every helix span. A helix span
`[helix-5p.start, helix-3p.end]` contains its inner loop span. That loop span
contains each child helix span. A hairpin span contains only its closing pair
and its unpaired segment.

Paint a span only when the color is meant to cover the nested content as
well. Paint outer spans first and inner spans after them, so a child replaces
the parent color on the bases it owns.

The closing pair is the boundary of the inner loop and is also two bases of
the helix. In the segment partition those two bases belong to the helix
strands. The loop span includes them so the span is a closed interval; the
segment list is what assigns each base a single color.

```
exterior span
└─ helix span
   └─ loop span          (hairpin, bulge, interior, or multibranch)
      ├─ unpaired segments of this loop
      └─ child helix span
         └─ …
```

### Worked example

```
.((..((...))..((...)).))
```

Indices `0 .. 23`. Pairs `(1,23), (2,22), (5,11), (6,10), (14,20), (15,19)`.

| loop_id | kind | start | end | parent | closing helix | unpaired segments |
|:--------|:-----|------:|----:|:-------|:--------------|:------------------|
| 0 | exterior | 0 | 23 | — | — | `[0,0]` |
| 1 | multibranch | 2 | 22 | 0 | H0 | `[3,4]`, `[12,13]`, `[21,21]` |
| 2 | hairpin | 6 | 10 | 1 | H1 | `[7,9]` |
| 3 | hairpin | 15 | 19 | 1 | H2 | `[16,18]` |

| helix_id | 5' strand | 3' strand | span | outer loop | inner loop |
|:---------|:----------|:----------|:-----|:-----------|:-----------|
| H0 | `[1,2]` | `[22,23]` | `[1,23]` | 0 | 1 |
| H1 | `[5,6]` | `[10,11]` | `[5,11]` | 1 | 2 |
| H2 | `[14,15]` | `[19,20]` | `[14,20]` | 1 | 3 |

Segments in order, a partition of `0 .. 23`:

| id | kind | start | end | link |
|---:|:-----|------:|----:|:-----|
| 0 | exterior | 0 | 0 | loop 0 |
| 1 | helix-5p | 1 | 2 | H0 |
| 2 | multibranch | 3 | 4 | loop 1 |
| 3 | helix-5p | 5 | 6 | H1 |
| 4 | hairpin | 7 | 9 | loop 2 |
| 5 | helix-3p | 10 | 11 | H1 |
| 6 | multibranch | 12 | 13 | loop 1 |
| 7 | helix-5p | 14 | 15 | H2 |
| 8 | hairpin | 16 | 18 | loop 3 |
| 9 | helix-3p | 19 | 20 | H2 |
| 10 | multibranch | 21 | 21 | loop 1 |
| 11 | helix-3p | 22 | 23 | H0 |

Containment of spans: H0 contains loop 1, which contains H1 and H2, which
contain loops 2 and 3. Segment rows do not contain each other.

With `circular` false the optional strings are

```
naview_graph: E1(2:M2+2+1(2:P3,2:P3))
rag_tree:     (*,@,@)
rag_dual:     (*-@-@-)
```

`E1` is the one exterior base. The multibranch gaps are 2, 2, and 1. Each
child helix has 2 base pairs and each hairpin has 3 unpaired bases. The RAG
tree keeps the end vertex, the junction, and both hairpins, so the center is
the junction. The dual is a 3-cycle of those stems, with a self-loop on each
hairpin stem.

With `circular` true the exterior is reclassified. It has one helix and one
unpaired base, so it is a hairpin whose unpaired base sits outside that
helix, and the strings are

```
naview_graph: P1(2:M2+2+1(2:P3,2:P3))
rag_tree:     (@,@)
rag_dual:     (@=@)
```

The closing hairpin has one unpaired base, so the RAG rules drop it and drop
the helix of length 2 that led to it. The two gaps that lay beside that
helix merge. What remains is a degree-2 loop between two hairpins, and the
dual has two parallel strands plus a self-loop on each stem. Base 0 is the
`hairpin` segment `[0,0]` of the closing loop. The other segment rows stay,
and the exterior row is gone.

An interior loop is the same shape with one child helix and two non-empty
gaps, as in `((..((...))..))`: outer helix `[0,1]`/`[13,14]`, inner helix
`[4,5]`/`[9,10]`, interior gaps `[2,3]` (`side = 5p`) and `[11,12]`
(`side = 3p`), hairpin `[6,8]`. A bulge is that pattern with one gap empty,
as in `((..((...))))`: the 5' gap is a `bulge` segment and the 3' side emits
no segment.

## Graph notation

Both optional outputs are single strings with no spaces. They describe the
loop tree. The pictures of those strings are a separate layout.

### NAView graph

`naview_graph` is the loop tree with loop kinds, gap sizes, and helix
lengths. The root is the exterior when `circular` is false and the closing
loop when `circular` is true. Children are in backbone order, starting at the
gap that contains base `0`.

```
tree   := node
node   := kind gaps children
kind   := 'E' | 'C' | 'P' | 'B' | 'I' | 'M'
gaps   := integer                  hairpin, circle, or a one-gap root
        | integer ('+' integer)+   every other loop
children := empty | '(' branch (',' branch)* ')'
branch := integer ':' node
```

`E` is exterior, `C` is circle, `P` is hairpin, `B` is bulge, `I` is
interior, `M` is multibranch. Each integer in `gaps` is the number of real
unpaired bases in that gap. The root lists every incident helix as a child,
so a root with `d` helices has `d` gaps and `d` children. Any other loop
omits the parent helix from the child list and still lists every gap, so a
non-root with `d` incident helices has `d` gaps and `d - 1` children. The
integer on a branch is the number of base pairs in that helix. A hairpin or
a circle with no child is `P3` or `C24`. A zero gap is written `0`.

The same dot-bracket produces a different string when `circular` changes,
because the closing loop changes kind and the 5' and 3' tails become one
gap.

### RAG tree

`rag_tree` is the tree of Gan et al., *Bioinformatics* 20:1285–1291 (2004),
built from the NAView tree by the following contractions. Every dot-bracket
pair counts as complementary, so only the length filters below apply. Edge
length is the number of base pairs. A vertex is a **junction** when it has
three or more incident helices. The exterior is an **end** vertex when
`circular` is false. The closing loop of a circular molecule is not an end
vertex.

1. Delete every helix of length 1 and fuse the two loops it joined. The
   fused vertex sums their unpaired bases, and it is a junction if the
   fused degree is at least 3. It is an end vertex if either side was.
2. Repeat until no vertex qualifies. A vertex qualifies when it is not a
   junction, not an end vertex, and its unpaired count is at most 1.
   - Degree 1: delete the vertex and its helix. At the neighbor, merge
     the two gaps that flanked that helix into one gap, summing their
     unpaired bases.
   - Degree 2: delete the vertex and replace its helices of lengths `a`
     and `b` by one helix of length `a + b`.
3. A circular molecule that loses every helix is one vertex. A linear
   molecule keeps the end vertex even when no helix remains.

Leaf marks record the only distinction the dual needs. A degree-1 hairpin
(unpaired count greater than 1) is `@`. A degree-1 end vertex is `*`. Every
surviving vertex of degree 2 or more is written with parentheses.

The string is centered so it does not depend on which loop was the drawing
root or on where a circular molecule was cut:

- One vertex: `*` or `@`.
- The center of a tree is the one vertex of minimum eccentricity, or the
  two adjacent vertices when they share that eccentricity.
- One center: root the tree there. A leaf child is `*` or `@`. Any other
  child is the parenthesized form of its own sorted children. Sort the
  children lexicographically and write `(c1,c2,...,ck)`.
- Two centers: write the rooted form of each side of the central edge,
  joined by `-`, with the lexicographically smaller side first. Two
  single-vertex sides give `*-*`, `*-@`, or `@-@`.

### RAG dual

`rag_dual` is the dual from the same paper, built from the contracted graph
and its gaps. The tree string alone does not record how many strands a loop
still has, so a bulge and an interior can share one `rag_tree` and differ
here.

A dual vertex is a helix that survived. A helix whose leaf side is `@` has
a self-loop. Each surviving loop adds one dual edge per non-empty gap, in
backbone order, between the consecutive helices of that loop. An empty gap
adds none. Two positive gaps between the same pair are two parallel edges.

```
dual  := leaf | '(' piece (sep piece)* sep? ')'
piece := dual
sep   := '-' | '=' | '.'
leaf  := '*' | '@'
```

`@` is a stem with one hairpin self-loop. `*` is a stem with none. `-` is
one strand, `=` is two parallel strands, and `.` is an empty gap. With three
or more pieces the separator before `)` closes the cycle and is required.
With two pieces the one separator is every strand between that pair.

| Tree | Dual | Reading |
|:-----|:-----|:--------|
| `*` | `*` | No helix. |
| `*-@` | `@` | One stem, one hairpin self-loop. |
| `@-@` | `@@` | One stem, a hairpin self-loop at each end. |
| `(*,@,@)` | `(*-@-@-)` | Three stems in a cycle, self-loops on the hairpin stems. |
| `(@,@)` | `(@=@)` | Two hairpin stems joined by two strands. |
| `(*,@)` and one empty gap | `(*-@)` | Two stems joined by the one strand of a bulge. |

## Graph layout

A graph layout is a second picture, independent of the nucleotide coordinates.
It draws `naview_graph`, `rag_tree`, or `rag_dual`. Scale, origin, rotation,
and which way is up are arbitrary, in the same sense as the base layout. The
reference frame used below has edge length 1, the root at the origin, `+x`
to the right, `+y` up, and angles increasing counterclockwise from `+x`.

Each vertex is a point `(x, y)`. Each drawn stroke is a line `(x0, y0)` to
`(x1, y1)`. A tree edge and a single dual strand are one line. A parallel
pair is two lines. A hairpin self-loop is three lines. Every stroke has two
endpoints.

### Tree lines

`naview_graph` and `rag_tree` use the same placer. `naview_graph` is rooted
at its notation root, and its children stay in backbone order. `rag_tree` is
rooted at the center already used to write the string, and its children stay
in that written order.

A single vertex sits at the origin and emits no line. A bicenter (`*-@`,
`@-@`, or `*-*`) places its two vertices at `(0, 0)` and `(1, 0)` and emits
that one line. Each side then places any further children in the open
half-turn that faces away from the other vertex: `(π/2, 3π/2)` at `(0, 0)`,
and `(-π/2, π/2)` at `(1, 0)`.

Otherwise the root is at the origin and owns the wedge `[0, 2π)`. A leaf
weighs 1. Any other vertex weighs the sum of its children. To place the
children of a vertex at `q` that owns `[α, β)`:

- Let `W` be that vertex's weight, and set `t = α`.
- For each child, in order, the wedge width is `(β - α) · weight(child) / W`.
  The child sits at distance 1 from `q` along the bisector `t + width/2`:

  ```
  child = q + (cos(t + width/2), sin(t + width/2))
  ```

- Emit the line from `q` to that child. The child owns `[t, t + width)` for
  its own children. Then advance `t` by `width`.

At the root this width is `2π · weight(child) / weight(root)`.

### Dual lines

Walk `rag_dual` in written order and collect the stem vertices (`*` or `@`).
Place `m` stems on the unit circle:

```
stem k = (cos(2π k / m), sin(2π k / m))
```

`m = 0` emits no point and no line. `m = 1` (`@` or `@@`) sits at the origin.

A separator `-` between two distinct stems is one line between their points.
A separator `=` is two lines, offset by `±0.2` along the unit perpendicular
of the chord. With `d` the unit vector from the earlier stem to the later
one, that perpendicular is `p = (-d_y, d_x)`, and the endpoints are the stem
points shifted by `±0.2 p`. A closing separator before `)` is the same kind
of line from the last piece back to the first. A `.` emits no line.

Each `@` is a triangle of three lines beside the stem. Let `u` be the
outward unit radial (the stem point itself when it lies on the unit circle,
and `(1, 0)` when the stem is the origin). Let `p = (-u_y, u_x)`. The
corners are

```
A = v + 0.45 u
B = v + 0.25 u + 0.25 p
C = v + 0.25 u - 0.25 p
```

and the lines are `A→B`, `B→C`, and `C→A`. The second loop of `@@` uses
`-u` in place of `u`.

### Lines for the worked example

Linear `rag_tree` `(*,@,@)`. The junction is the origin. The three leaves
weigh 1 each, so their bisectors are `π/3`, `π`, and `5π/3`, a third of a
turn apart.

| Vertex | `(x, y)` |
|:-------|:---------|
| junction | `(0, 0)` |
| `*` | `(1/2, √3/2)` |
| first `@` | `(-1, 0)` |
| second `@` | `(1/2, -√3/2)` |

| Line | Start | End |
|:-----|:------|:----|
| junction–`*` | `(0, 0)` | `(1/2, √3/2)` |
| junction–first `@` | `(0, 0)` | `(-1, 0)` |
| junction–second `@` | `(0, 0)` | `(1/2, -√3/2)` |

Linear `rag_dual` `(*-@-@-)`. Three stems on the unit circle, three chords,
and a triangle on each `@`.

| Stem | Mark | `(x, y)` |
|:-----|:-----|:---------|
| 0 | `*` | `(1, 0)` |
| 1 | `@` | `(-1/2, √3/2)` |
| 2 | `@` | `(-1/2, -√3/2)` |

| Line | Start | End |
|:-----|:------|:----|
| `0-1` | `(1, 0)` | `(-1/2, √3/2)` |
| `1-2` | `(-1/2, √3/2)` | `(-1/2, -√3/2)` |
| `2-0` | `(-1/2, -√3/2)` | `(1, 0)` |
| `@` on stem 1 | `(-0.725, 0.725√3)` | `((-5-√3)/8, (5√3-1)/8)` |
| `@` on stem 1 | `((-5-√3)/8, (5√3-1)/8)` | `((-5+√3)/8, (5√3+1)/8)` |
| `@` on stem 1 | `((-5+√3)/8, (5√3+1)/8)` | `(-0.725, 0.725√3)` |
| `@` on stem 2 | `(-0.725, -0.725√3)` | `((-5+√3)/8, (-5√3-1)/8)` |
| `@` on stem 2 | `((-5+√3)/8, (-5√3-1)/8)` | `((-5-√3)/8, (-5√3+1)/8)` |
| `@` on stem 2 | `((-5-√3)/8, (-5√3+1)/8)` | `(-0.725, -0.725√3)` |

Circular `rag_tree` `(@,@)`. The degree-2 vertex is the origin. Its two
hairpins sit at `(0, 1)` and `(0, -1)`. The lines are `(0, 0)–(0, 1)` and
`(0, 0)–(0, -1)`.

Circular `rag_dual` `(@=@)`. The stems are `(1, 0)` and `(-1, 0)`. The
chord direction is `(-1, 0)`, so `p = (0, -1)`, and the two strands are

| Line | Start | End |
|:-----|:------|:----|
| `+0.2 p` | `(1, -0.2)` | `(-1, -0.2)` |
| `-0.2 p` | `(1, 0.2)` | `(-1, 0.2)` |

Each `@` adds its triangle. On `(1, 0)` the corners are `(1.45, 0)`,
`(1.25, 0.25)`, and `(1.25, -0.25)`. On `(-1, 0)` they are `(-1.45, 0)`,
`(-1.25, -0.25)`, and `(-1.25, 0.25)`. The six triangle sides run corner to
corner in that order, returning to the first corner.

`naview_graph` uses this same tree placer, rooted at `E` or at `P1`, with
wedges in backbone order. For `E1(2:M2+2+1(2:P3,2:P3))` the exterior is the
origin, the multibranch is at `(-1, 0)`, and the hairpins are at `(-1, 1)`
and `(-1, -1)`. The lines are `(0, 0)–(-1, 0)`, `(-1, 0)–(-1, 1)`, and
`(-1, 0)–(-1, -1)`. `P1(2:M2+2+1(2:P3,2:P3))` has the same three lines,
because that root also has one child and the child has two equal leaves.

## Algorithm

### Pairs

Scan left to right with a stack of open-bracket indices.

- `(` pushes `i`.
- `)` pops `j` and records the pair both ways. Reject the string if the stack
  is empty.
- `.` records no partner.

Reject the string if the stack is not empty at the end. Partner lookup is
then: `mate[i] = j` or `mate[i]` absent.

### Helices

A helix is a maximal run of stacked pairs. Walk `i` from `0` to `n-1`. When
`i` is paired with `j > i` and the pair is not yet in a helix:

- `start1 = i`, `end2 = j`
- While `i + 1 < j - 1` and base `i + 1` is paired with `j - 1`, set
  `i = i + 1` and `j = j - 1`.
- `end1 = i` and `start2 = j`. The test `i + 1 < j - 1` keeps the two
  strands from crossing, including a helix that ends in a pair with no
  unpaired base between them.

The 5' strand is `[start1, end1]`. The 3' strand is `[start2, end2]`. Pairing
inside the helix is `start1+k` with `end2-k`.

### Loop tree

Loops are the nodes of a tree. Helices are the edges. Build the tree by a
circular walk that starts on the exterior, which is the closing loop when
`circular` is true. When `circular` is false, a phantom slot that is not a
base sits between base `n-1` and base `0`, and the walk crosses it. When
`circular` is true, the walk crosses a real backbone bond and there is no
phantom base.

At each loop, walk forward from its entry base. On reaching the outer end of
a helix that has not yet been entered, mark that helix entered, build the
loop on its inner side by recursion, and link the two loops through the
helix. Then jump to the partner and continue, so the walk does not travel
down the helix base by base. The inner loop is entered at the first base
after the inner end of the strand just used (`end1+1`, or `end2+1` when the
walk hit the 3' strand first). Wrapping past `n-1` continues at `0`.

The link stores, on each of the two loops, the pair that faces that loop:

- On the outer loop, the facing pair is the outer pair of the helix.
- On the inner loop, the facing pair is the inner pair of the helix. That
  pair is the inner loop's `(start, end)`.

Discovery order is preorder: loop `0` is created first, and the helices
around a loop are entered in backbone order. That order assigns `loop_id` and
is the tie break for the root below. After the tree exists, reclassify loop
`0` when `circular` is true, using the kind table and the real closing bond.
Its 5' tail and 3' tail are one gap.

### Root

The drawing root is the most central highly branched loop.

The depth of a loop is the number of helices on the shortest path from that
loop to a leaf. A leaf is a loop with one helix (every hairpin, including a closing loop
that has only one helix). An all-unpaired molecule is handled in
[No pairs](#no-pairs) and does not use this rule.

Choose the loop with the largest number of incident helices. If several loops
share that number, choose the one with the largest depth. If several still
remain, keep the one with the smallest `loop_id`.

### Circle diagram

When `circular` is false, reserve `n+1` equal angles, one per base plus the
phantom slot between the 3' end and the 5' end. The phantom has angle `0`.
Base `i` has angle

```
θ(i) = 2π * (i + 1) / (n + 1)
```

When `circular` is true there is no phantom. Base `i` has angle

```
θ(i) = 2π * i / n
```

with `θ(0) = 0` for `n = 1`. Neighboring bases, including base `n-1` and
base `0`, are one step apart on this diagram. Stem directions therefore
differ from the linear layout, and the coordinates differ by more than a
scale, rotation, or translation.

Stem directions are the normals of the chords of this diagram. They depend on
sequence position, so the same helix gets the same direction in a related
structure. For a pair that faces the current loop, with bases `p` and `q`:

```
(xs, ys) = (-sin θ(p), cos θ(p))
(xe, ye) = (-sin θ(q), cos θ(q))
radial   = normalize(ye - ys, xs - xe)
stem_angle = atan2(radial_y, radial_x)     in [0, 2π)
```

`radial` points from the loop center out along the helix.

The step count from a facing base `e` forward to the next facing base `s`
is the number of backbone links in that gap. Let `period` be `n + 1` when
`circular` is false and `n` when `circular` is true.

```
steps(e, s) = (s - e)            if s > e
            = s + period - e     otherwise
```

When `circular` is false the phantom counts as one link, so a gap that
crosses the join has `steps - 2` real unpaired bases, and any other gap has
`steps - 1`. When `circular` is true every link is a real backbone bond, and
every gap has `steps - 1` real unpaired bases. The closing bond alone is
`steps(n - 1, 0) = 1`.

### Radius

Backbone links and base pairs have target length 1. Around a loop, let each
gap between consecutive facing pairs (in backbone order, wrapping from the
last back to the first) have circle-diagram angle `dt > 0` and link count
`ci`.

- A gap that is not extruded has `ci = steps(e, s)`.
- An extruded gap has `ci = 2` when `dt ≤ π/2`, and `ci = 1.5` otherwise.
  Those constants stand in for the extruded arc, so that gap does not force
  the circle to grow.

The radius is

```
r = Σ dt * (1/ci + 1)  /  Σ (dt² / ci)
```

summed over the gaps of that loop. If `r < √2 / 2`, use `√2 / 2`.

Let `lencut = 0.5`. Among gaps that are not extruded and have `ci > 1`, let
`mindit` be the smallest `dt/ci`. If `r * mindit < lencut`, mark that gap
extruded and compute `r` again. Stop when every remaining eligible gap has
`r * dt/ci ≥ lencut`, or when no eligible gap remains. A gap with `ci = 1`
is a backbone bond between two helices and is never extruded by this test.

### Placing a loop

Place the root first, then each child after its helix has been drawn.

1. Compute `radial` and `stem_angle` for every helix that faces the loop.
2. Compute `r`, extruding gaps as above.
3. Set the center. The root center is `(0, 0)`. Any other center is the
   midpoint of the anchor pair (already placed by the parent) minus
   `r * radial` of that pair, so the anchor pair lies on this loop's circle
   and the helix points outward.
4. Place every facing pair. Then place the unpaired bases. Then, for each
   child helix, draw the helix and place the child loop.

Two facing pairs are **chained** when the gap between them is extruded or
`steps(e, s) = 1`. Chained pairs are laid out as one rigid run. Split the
helices around the loop into maximal chained runs. If the whole loop is one
chained run, break the run at the largest circle-diagram gap so the run has
two ends.

For each run:

- The middle of the run is the anchor helix when the run contains the helix
  this loop was entered through. Otherwise, with the helices of the run
  numbered `0 .. count-1` in backbone order, the middle is index
  `floor((count - 1) / 2)`.
- The anchor pair keeps the coordinates written when its helix was drawn.
  Every other middle pair is placed on the circle, centered on its
  `stem_angle`, with chord length 1. With `α = asin(1 / (2r))`:

  ```
  start base = center + r * (cos(stem_angle - α), sin(stem_angle - α))
  end base   = center + r * (cos(stem_angle + α), sin(stem_angle + α))
  ```

- Walk outward to each end of the run. The next pair is offset from the
  previous pair by length `rl` along the unit perpendicular to the bisector
  of the two stem angles (the perpendicular that continues away from the
  middle). Its own two bases are then placed 1 apart, perpendicular to its
  `radial`. Use `rl = 1` when the gap just crossed is not extruded. When it
  is extruded, use `rl = 2` if that gap's `dt ≤ π/2`, and `rl = 1.5`
  otherwise.
- A run of one pair is already on the circle, and the single run that goes
  all the way around the loop is left as built. Every other run is moved so
  the midpoint of its two extreme bases lies on the circle, on the side of
  the shorter arc between those bases. When the run contains the anchor
  pair, move the center and leave the bases of the run in place. Otherwise
  move those bases and leave the center in place.

After every facing pair has coordinates, compare the counterclockwise
geometric angle from each end base to the next start base with that gap's
circle-diagram angle. If they differ by more than `π`, the gap has crossed
itself. When that gap is not yet extruded and it contains at least one real base,
mark it extruded and return to step 2. In the linear layout the historical
test is `s ≠ e + 1`, which also matches a phantom-only join; keep that test
when `circular` is false. When `circular` is true, base `0` is the next base
after base `n-1`, and a closing bond with no unpaired base is not extruded
by this test.

Unpaired bases in a gap that is not extruded are spaced at equal angle along
that counterclockwise arc. The radius varies linearly with the angle between
the two bases' distances from the center. When `circular` is false, a gap
that crosses the phantom spends one angular station on the phantom and emits
no base for it.

An extruded gap is drawn as a circular arc of backbone length 1 when that
arc fits. Let the chord between the current endpoints have length `b`, and
let `m` be the number of backbone links still to be drawn (`m ≥ 1`). If
`b ≥ m`, place those bases evenly on the straight chord. Otherwise find a
center at signed distance `h` from the chord midpoint such that the remaining
sides have length 1 and the polygon closes:

```
ρ     = sqrt(h² + (b/2)²)
θside = acos(1 - 0.5 / ρ²)
φ     = acos(h / ρ)
error = θside * m + 2φ - 2π
```

Bisect `h`. The upper bound starts at `m / π`. The lower bound starts at
`0` when `b < 1`, and otherwise at `-(m / π) - b / (m + 1e-6 - b)`. Positive
`h` puts the center on the clockwise-normal side of the directed chord;
negative `h` puts it on the other side. Stop at `|error| ≤ 1e-4` or after
500 splits. If a split is undefined (`|1 - 0.5/ρ²| > 1`) or the iteration
does not converge, place the bases evenly on the straight chord.

If the arc leaves either stem in a direction that turns more than `π` from
that stem's circle-diagram angle, the arc does not fit. Move one base off
each end onto a unit step at an angle halfway from that stem angle toward
the gap's angle bisector (and at least `0.5` rad from the stem angle),
shorten the unsolved run by those two bases, and solve the arc again. Repeat
until the arc fits or a single link remains. A gap of exactly one
intermediate base is just the arc between its two stems, with no unit-step
peeling.

### Helix coordinates

The facing pair of a helix is already placed. Number the 5' strand
`start, start+1, …`. Base `start+k` and its partner are placed at

```
facing_5p + k * radial
facing_3p + k * radial
```

with the `radial` of the outer loop, which points out along the helix. Then
place the inner loop, anchored on the inner pair.

### No pairs

If the string contains no brackets, the output is one loop of kind
`exterior` when `circular` is false and kind `circle` when `circular` is
true, with span `[0, n-1]`, one segment of that same kind covering the span,
and no helix. The strings are `E<n>` or `C<n>`, and both RAG strings are
`*`.

Coordinates are a regular layout of bond length 1. One base sits at the
origin. For `n ≥ 2` the bases sit on a circle at equal angles with neighbor
chord 1. An unpaired chain uses this circle for either value of `circular`,
because it has no stem whose direction a phantom could move:

```
r = 1 / (2 * sin(π / n))
(x, y) = r * (cos(2π i / n), sin(2π i / n))
```

The historical program invented a fictitious pair between the first and last
base so that a pair always existed. This specification draws the circle
instead.

## Coloring

Segments are the ranges to color in a post-processing step.

- Color a base from its single segment. No segment contains another, so the
  order of segments does not matter.
- Color a stem by both of its strand segments.
- Color only the unpaired bases of a loop by the segments that carry that
  `loop_id`.
- Color a loop and everything nested in its closing pair by that pair's
  span, then repaint every properly contained loop or helix span that should
  keep a different color. Children are painted after parents. The closing
  loop of a circular molecule has no closing pair; color it from its
  segments. A segment with `start > end` is one of those regions: it covers
  `start .. n-1` and `0 .. end`.

Loop kind and segment kind use the same names, so a segment can be colored
from its own `kind` without looking up the loop. `helix-5p` and `helix-3p`
are the two strand colors of a stem; they share `helix_id` when the whole
stem is one color.
