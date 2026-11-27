'''
Write python function to parse a circular dot-bracket string with (.) only.
Output is dictionary with the following features:

# Features of the complete circle
nr_bases
nr_paired_bases # number of parentheses '(' or ')'
nr_unpaired_bases # number of dots '.'
pct_paired_bases # 100*nr_paired_bases/nr_bases

# Features of the longest stem
longest_stem_total_nr_bases # including unpaired
longest_stem_nr_paired_bases # counting only paired
longest_stem_pct_paired_bases # 100*longest_stem_nr_paired_bases/longest_stem_total_nr_bases
longest_stem_pct_total_bases # 100*longest_stem_total_nr_bases/nr_bases

# Features of the longest composite stem
longest_composite_stem_total_nr_bases # including unpaired
longest_composite_stem_nr_paired_bases
longest_composite_stem_pct_paired_bases # 100*longest_composite_stem_nr_paired_bases/longest_composite_stem_total_nr_bases
longest_composite_stem_pct_total_bases # 100*longest_composite_stem_total_nr_bases/nr_bases

# Bulginess
max_bulge_len # Max bulge (defined as consecutive dots)
nr_excess_bulges # Number of bulges > training set median
total_nr_bulge_bases # Number of bases in bulges
pct_bulge # 100*total_nr_bulge_bases/nr_bases
excess_bulge_nr_bases # Number of bases in excess bulges
pct_excess_bulges # 100*excess_bulge_nr_bases/nr_bases
bulge_to_stem_ratio # total_nr_bulge_bases/nr_paired_bases
bulge_length_L2_norm # sqrt(sum of bulge length squared)
nr_bases_true_bulges # unpaired on one strand only
nr_bases_symmetric_loops # unpaired on both strands
pct_true_bulges # 100*nr_bases_true_bulges/nr_bases
pct_symmetric_loops # 100*nr_bases_symmetric_loops/nr_bases

# ROD v2 structure. pS equals pct_paired_bases.
# pH, pI, pB, pF are hairpin, interior, bulge, and bifurcation loops.
# With no pairs those five, pRod, nBranch, pLongestRod and the spikiness
# percentages are 0.
pS
pH
pI
pB
pF
pRod
nBranch
pLongestRod
pOffRod # 100 * pairs off the best hairpin-to-hairpin path / nr_paired_bases
pOnRod # 100 - pOffRod; 0 when nothing is paired
pSplitStem # 100 * pairs outside the longest composite stem / nr_paired_bases
pLongestPairShare # 100 - pSplitStem; 0 when nothing is paired
pJunction # 100 * non-heaviest spoke pairs / spoke pairs, weighted over bifurcations
pDominantSpoke # 100 - pJunction; 0 when nothing is paired
bKeto # abs(GT-AC)/(GT+AC) of paired bases; 0 with no pairs or no sequence
uKeto # same skew of unpaired bases; 0 with no sequence
'''

import math
from typing import Dict, List, Optional, Tuple


class _Loop:
	def __init__(self, kind: str, unpaired: int, incident: List[Tuple[int, int]]):
		self.kind = kind
		self.unpaired = unpaired
		self.incident = incident
		self.anchor = kind in ("H", "F")


def _pair_tree(
	dot_bracket: str, pair_map: Dict[int, int]
) -> Tuple[
	Dict[Tuple[int, int], Optional[Tuple[int, int]]],
	Dict[Tuple[int, int], List[Tuple[int, int]]],
	List[Tuple[int, int]],
]:
	parent: Dict[Tuple[int, int], Optional[Tuple[int, int]]] = {}
	opening_stack: List[int] = []
	for index, char in enumerate(dot_bracket):
		if char == "(":
			opening_stack.append(index)
		elif char == ")":
			opening = opening_stack.pop()
			pair = (opening, index)
			if opening_stack:
				parent_opening = opening_stack[-1]
				parent[pair] = (parent_opening, pair_map[parent_opening])
			else:
				parent[pair] = None
	children: Dict[Tuple[int, int], List[Tuple[int, int]]] = {pair: [] for pair in parent}
	for pair, enclosing in parent.items():
		if enclosing is not None:
			children[enclosing].append(pair)
	for kids in children.values():
		kids.sort()
	top_level = sorted(pair for pair, enclosing in parent.items() if enclosing is None)
	return parent, children, top_level


def _arcs_inside(opening: int, closing: int, kids: List[Tuple[int, int]]) -> List[int]:
	arcs: List[int] = []
	cursor = opening + 1
	for child_open, child_close in kids:
		arcs.append(child_open - cursor)
		cursor = child_close + 1
	arcs.append(closing - cursor)
	return arcs


def _exterior_arcs(top_level: List[Tuple[int, int]], length: int) -> List[int]:
	arcs: List[int] = []
	count = len(top_level)
	for index, (_opening, closing) in enumerate(top_level):
		next_opening = top_level[(index + 1) % count][0]
		if index + 1 < count:
			arcs.append(next_opening - closing - 1)
		else:
			arcs.append((length - closing - 1) + next_opening)
	return arcs


def _loop_kind(arcs: List[int], degree: int) -> Tuple[str, int]:
	unpaired = sum(arcs)
	if degree == 1:
		return "H", unpaired
	if degree == 2:
		occupied = sum(1 for arc in arcs if arc > 0)
		if occupied == 2:
			return "I", unpaired
		if occupied == 1:
			return "B", unpaired
		return "stack", 0
	return "F", unpaired


def _rod_percentages(dot_bracket: str, pair_map: Dict[int, int]) -> Dict[str, float]:
	"""Circular loop states and the rod backbone. All zero when nothing is paired."""
	length = len(dot_bracket)
	zeros = {
		"pS": 0.0,
		"pH": 0.0,
		"pI": 0.0,
		"pB": 0.0,
		"pF": 0.0,
		"pRod": 0.0,
		"nBranch": 0.0,
		"pLongestRod": 0.0,
		"pOffRod": 0.0,
		"pOnRod": 0.0,
		"pJunction": 0.0,
		"pDominantSpoke": 0.0,
	}
	if length == 0 or not pair_map:
		return zeros

	_parent, children, top_level = _pair_tree(dot_bracket, pair_map)
	inside: Dict[Tuple[int, int], _Loop] = {}
	outside: Dict[Tuple[int, int], _Loop] = {}
	loops: List[_Loop] = []

	for pair, kids in children.items():
		opening, closing = pair
		kind, unpaired = _loop_kind(_arcs_inside(opening, closing, kids), 1 + len(kids))
		loop = _Loop(kind, unpaired, [pair] + kids)
		loops.append(loop)
		inside[pair] = loop

	exterior_kind, exterior_unpaired = _loop_kind(_exterior_arcs(top_level, length), len(top_level))
	exterior = _Loop(exterior_kind, exterior_unpaired, list(top_level))
	loops.append(exterior)
	for pair in top_level:
		outside[pair] = exterior
	for pair, enclosing in _parent.items():
		if enclosing is not None:
			outside[pair] = inside[enclosing]

	hairpin_bases = 0
	interior_bases = 0
	bulge_bases = 0
	bifurcation_bases = 0
	n_branch = 0
	for loop in loops:
		if loop.kind == "H":
			hairpin_bases += loop.unpaired
		elif loop.kind == "I":
			interior_bases += loop.unpaired
		elif loop.kind == "B":
			bulge_bases += loop.unpaired
		elif loop.kind == "F":
			bifurcation_bases += loop.unpaired
			n_branch += 1

	def other_side(loop: _Loop, pair: Tuple[int, int]) -> _Loop:
		if inside[pair] is loop:
			return outside[pair]
		return inside[pair]

	# A stretch joins two hairpin or bifurcation loops through bulges and interior loops.
	used = set()
	adjacent: Dict[int, List[Tuple[_Loop, int, int]]] = {}
	anchors = [loop for loop in loops if loop.anchor]
	for anchor in anchors:
		adjacent[id(anchor)] = []
	for anchor in anchors:
		for pair in anchor.incident:
			if pair in used:
				continue
			collected = [pair]
			used.add(pair)
			unpaired = 0
			arrived = pair
			current = other_side(anchor, pair)
			while not current.anchor:
				if len(current.incident) != 2:
					raise ValueError("degree-2 loop does not have two pairs")
				unpaired += current.unpaired
				first, second = current.incident
				nxt = second if arrived == first else first
				collected.append(nxt)
				used.add(nxt)
				arrived = nxt
				current = other_side(current, nxt)
			paired_bases = 2 * len(collected)
			total_bases = paired_bases + unpaired
			adjacent[id(anchor)].append((current, paired_bases, total_bases))
			adjacent[id(current)].append((anchor, paired_bases, total_bases))

	if len(used) != len(children):
		raise ValueError("a base pair is missing from the rod backbone")

	hairpins = [loop for loop in anchors if loop.kind == "H"]
	if len(hairpins) < 2:
		raise ValueError("a paired structure has fewer than two hairpin loops")

	def path_edges(source: _Loop, target: _Loop) -> List[Tuple[int, int]]:
		parent: Dict[int, Optional[_Loop]] = {id(source): None}
		edge_to: Dict[int, Tuple[int, int]] = {}
		stack = [source]
		while stack:
			node = stack.pop()
			if node is target:
				break
			for other, paired_bases, total_bases in adjacent[id(node)]:
				if id(other) in parent:
					continue
				parent[id(other)] = node
				edge_to[id(other)] = (paired_bases, total_bases)
				stack.append(other)
		edges: List[Tuple[int, int]] = []
		node = target
		while node is not source:
			edges.append(edge_to[id(node)])
			parent_loop = parent[id(node)]
			if parent_loop is None:
				raise ValueError("hairpin loops are not connected")
			node = parent_loop
		return edges

	best_key: Optional[Tuple[int, int, int]] = None
	best_paired = 0
	best_longest_pairs = 0
	for left in range(len(hairpins)):
		for right in range(left + 1, len(hairpins)):
			edges = path_edges(hairpins[left], hairpins[right])
			paired = sum(item[0] for item in edges)
			longest_total, longest_pairs = max((item[1], item[0]) for item in edges)
			key = (paired, longest_total, longest_pairs)
			if best_key is None or key > best_key:
				best_key = key
				best_paired = paired
				best_longest_pairs = longest_pairs

	def percent(count: int) -> float:
		return 100.0 * count / length

	paired_count = len(pair_map)
	junction_weight = 0
	junction_discarded = 0
	for loop in loops:
		if loop.kind != "F":
			continue
		spokes = [paired for _other, paired, _total in adjacent[id(loop)]]
		weight = sum(spokes)
		if weight <= 0:
			continue
		junction_weight += weight
		junction_discarded += weight - max(spokes)
	if paired_count == 0:
		off_rod = 0.0
	else:
		off_rod = 100.0 * (paired_count - best_paired) / paired_count
	if junction_weight <= 0:
		junction = 0.0
	else:
		junction = 100.0 * junction_discarded / junction_weight

	return {
		"pS": percent(paired_count),
		"pH": percent(hairpin_bases),
		"pI": percent(interior_bases),
		"pB": percent(bulge_bases),
		"pF": percent(bifurcation_bases),
		"pRod": percent(best_paired),
		"nBranch": float(n_branch),
		"pLongestRod": percent(best_longest_pairs),
		"pOffRod": off_rod,
		"pOnRod": 0.0 if paired_count == 0 else 100.0 - off_rod,
		"pJunction": junction,
		"pDominantSpoke": 0.0 if paired_count == 0 else 100.0 - junction,
	}


def _keto_symbol(base: str) -> Optional[str]:
	char = base.upper()
	if char == "U":
		char = "T"
	if char in ("G", "T"):
		return "keto"
	if char in ("A", "C"):
		return "amino"
	return None


def _absolute_skew(symbols: List[Optional[str]]) -> float:
	gt = sum(1 for symbol in symbols if symbol == "keto")
	ac = sum(1 for symbol in symbols if symbol == "amino")
	if gt + ac == 0:
		return 0.0
	return abs(gt - ac) / (gt + ac)


def _keto_skews(dot_bracket: str, sequence: str) -> Dict[str, float]:
	if sequence == "":
		return {"bKeto": 0.0, "uKeto": 0.0}
	if len(sequence) != len(dot_bracket):
		raise ValueError(
			"sequence length %d does not match dot-bracket length %d"
			% (len(sequence), len(dot_bracket))
		)
	paired: List[Optional[str]] = []
	unpaired: List[Optional[str]] = []
	for index, char in enumerate(dot_bracket):
		symbol = _keto_symbol(sequence[index])
		if char == ".":
			unpaired.append(symbol)
		else:
			paired.append(symbol)
	bound = 0.0 if not paired else _absolute_skew(paired)
	return {"bKeto": bound, "uKeto": _absolute_skew(unpaired)}


def extract_circular_rna_features(
	dot_bracket: str,
	max_unpaired_gap: int = 5,
	train_median_bulge_len: float = 1.0,
	sequence: str = "",
) -> Dict[str, float]:
	"""
	Parses a 3-symbol dot-bracket string ('(', ')', '.') representing a circular RNA
	and extracts RCE stem features plus the ROD v2 loop, backbone, and keto features.

	Args:
		dot_bracket: Circular secondary structure in 3-symbol dot-bracket notation.
		max_unpaired_gap: Max consecutive unpaired nucleotides allowed per strand
			when bridging helices into composite stems and detecting bulges/loops.
		train_median_bulge_len: Median bulge run length from the training set
			used to define excess bulges (> train_median_bulge_len).
		sequence: Nucleotide string aligned with dot_bracket. Empty leaves bKeto
			and uKeto at 0. U counts as T. Letters other than A, C, G, T, U are skipped.

	Returns:
		Dictionary containing the RCE features and the ROD v2 structure features.
	"""
	n = len(dot_bracket)

	# Zero-initialized output dictionary for clean edge-case handling
	features: Dict[str, float] = {
		# Features of the complete circle
		"nr_bases": float(n),
		"nr_paired_bases": 0.0,
		"nr_unpaired_bases": 0.0,
		"pct_paired_bases": 0.0,
		# Features of the longest strict stem (no internal unpaired bases)
		"longest_stem_total_nr_bases": 0.0,
		"longest_stem_nr_paired_bases": 0.0,
		"longest_stem_pct_paired_bases": 0.0,
		"longest_stem_pct_total_bases": 0.0,
		# Features of the longest composite stem (bridges bulges/internal loops)
		"longest_composite_stem_total_nr_bases": 0.0,
		"longest_composite_stem_nr_paired_bases": 0.0,
		"longest_composite_stem_pct_paired_bases": 0.0,
		"longest_composite_stem_pct_total_bases": 0.0,
		# Bulginess
		"max_bulge_len": 0.0,
		"nr_excess_bulges": 0.0,
		"total_nr_bulge_bases": 0.0,
		"pct_bulge": 0.0,
		"excess_bulge_nr_bases": 0.0,
		"pct_excess_bulges": 0.0,
		"bulge_to_stem_ratio": 0.0,
		"bulge_length_L2_norm": 0.0,
		"nr_bases_true_bulges": 0.0,
		"nr_bases_symmetric_loops": 0.0,
		"pct_true_bulges": 0.0,
		"pct_symmetric_loops": 0.0,
		# ROD v2 structure. pS matches pct_paired_bases.
		"pS": 0.0,
		"pH": 0.0,
		"pI": 0.0,
		"pB": 0.0,
		"pF": 0.0,
		"pRod": 0.0,
		"nBranch": 0.0,
		"pLongestRod": 0.0,
		"pOffRod": 0.0,
		"pOnRod": 0.0,
		"pSplitStem": 0.0,
		"pLongestPairShare": 0.0,
		"pJunction": 0.0,
		"pDominantSpoke": 0.0,
		"bKeto": 0.0,
		"uKeto": 0.0,
	}

	if n == 0:
		return features

	# 1. Stack-based pairing extraction
	stack: List[int] = []
	pair_map: Dict[int, int] = {}

	for i, char in enumerate(dot_bracket):
		if char == "(":
			stack.append(i)
		elif char == ")":
			if not stack:
				raise ValueError(f"Unmatched ')' at index {i}")
			j = stack.pop()
			pair_map[j] = i
			pair_map[i] = j
		elif char != ".":
			raise ValueError(f"Invalid character '{char}' at index {i}")

	if stack:
		raise ValueError(f"Unmatched '(' at indices {stack}")

	nr_paired = len(pair_map)
	nr_unpaired = n - nr_paired

	features["nr_paired_bases"] = float(nr_paired)
	features["nr_unpaired_bases"] = float(nr_unpaired)
	features["pct_paired_bases"] = 100.0 * nr_paired / n
	features.update(_rod_percentages(dot_bracket, pair_map))
	features.update(_keto_skews(dot_bracket, sequence))

	canonical_pairs: List[Tuple[int, int]] = sorted(
		(i, j) for i, j in pair_map.items() if i < j
	)
	if not canonical_pairs:
		return features

	# Helper: collect circular indices strictly between start and end moving +1
	def circular_between(start: int, end: int) -> List[int]:
		indices = []
		curr = (start + 1) % n
		while curr != end:
			indices.append(curr)
			curr = (curr + 1) % n
		return indices

	# Helper: find next inward base pair bridging up to `gap` unpaired dots per strand
	def next_inward_step(
		i: int, j: int, gap: int
	) -> Optional[Tuple[Tuple[int, int], List[int], List[int]]]:
		for step_5 in range(1, gap + 2):
			ni = (i + step_5) % n
			if dot_bracket[ni] == ".":
				continue
			for step_3 in range(1, gap + 2):
				nj = (j - step_3) % n
				if dot_bracket[nj] == ".":
					continue
				if pair_map.get(ni) == nj and ni != j and nj != i:
					unpaired_5 = circular_between(i, ni)
					unpaired_3 = circular_between(nj, j)
					return (ni, nj), unpaired_5, unpaired_3
				return None
			return None
		return None

	# Helper: find previous outward base pair bridging up to `gap` unpaired dots per strand
	def prev_outward_pair(i: int, j: int, gap: int) -> Optional[Tuple[int, int]]:
		for step_5 in range(1, gap + 2):
			pi = (i - step_5) % n
			if dot_bracket[pi] == ".":
				continue
			for step_3 in range(1, gap + 2):
				pj = (j + step_3) % n
				if dot_bracket[pj] == ".":
					continue
				if pair_map.get(pi) == pj and pi != j and pj != i:
					return (pi, pj)
				return None
			return None
		return None

	# Helper: decompose circular structure into stems for a given gap tolerance
	def find_stems(gap: int) -> List[Dict]:
		visited = set()
		stems: List[Dict] = []

		for i, j in canonical_pairs:
			if (i, j) in visited or (j, i) in visited:
				continue

			# Walk outward to the outermost base pair of this stem
			curr = (i, j)
			while True:
				prev_pair = prev_outward_pair(*curr, gap=gap)
				if prev_pair and prev_pair not in visited and prev_pair != (i, j):
					curr = prev_pair
				else:
					break

			# Walk inward to collect pairs and loop segments
			stem_pairs: List[Tuple[int, int]] = []
			loop_segments: List[Tuple[List[int], List[int]]] = []

			while curr and curr not in visited:
				visited.add(curr)
				visited.add((curr[1], curr[0]))
				stem_pairs.append(curr)

				step_res = next_inward_step(*curr, gap=gap)
				if step_res:
					next_pair, u5, u3 = step_res
					if u5 or u3:
						loop_segments.append((u5, u3))
					curr = next_pair
				else:
					curr = None

			stems.append({
				"pairs": stem_pairs,
				"loop_segments": loop_segments,
			})

		return stems

	# 2. Longest strict stem (gap = 0 -> no internal unpaired bases)
	strict_stems = find_stems(gap=0)
	longest_strict = max(strict_stems, key=lambda s: len(s["pairs"]))
	strict_paired_bases = 2 * len(longest_strict["pairs"])
	strict_total_bases = strict_paired_bases  # always 0 unpaired in a strict stem

	features["longest_stem_total_nr_bases"] = float(strict_total_bases)
	features["longest_stem_nr_paired_bases"] = float(strict_paired_bases)
	features["longest_stem_pct_paired_bases"] = (
		100.0 * strict_paired_bases / strict_total_bases if strict_total_bases > 0 else 0.0
	)
	features["longest_stem_pct_total_bases"] = 100.0 * strict_total_bases / n

	# 3. Longest composite stem (gap = max_unpaired_gap)
	composite_stems = find_stems(gap=max_unpaired_gap)
	longest_comp = max(composite_stems, key=lambda s: len(s["pairs"]))
	comp_paired_bases = 2 * len(longest_comp["pairs"])
	comp_unpaired_bases = sum(
		len(u5) + len(u3) for u5, u3 in longest_comp["loop_segments"]
	)
	comp_total_bases = comp_paired_bases + comp_unpaired_bases

	features["longest_composite_stem_total_nr_bases"] = float(comp_total_bases)
	features["longest_composite_stem_nr_paired_bases"] = float(comp_paired_bases)
	features["longest_composite_stem_pct_paired_bases"] = (
		100.0 * comp_paired_bases / comp_total_bases if comp_total_bases > 0 else 0.0
	)
	features["longest_composite_stem_pct_total_bases"] = 100.0 * comp_total_bases / n
	if nr_paired > 0:
		features["pSplitStem"] = 100.0 * (1.0 - comp_paired_bases / nr_paired)
		features["pLongestPairShare"] = 100.0 - features["pSplitStem"]

	# 4. Bulginess across all composite stems in the circle
	bulge_runs: List[int] = []
	nr_true_bulge_bases = 0
	nr_symmetric_loop_bases = 0

	for stem in composite_stems:
		for u5, u3 in stem["loop_segments"]:
			len5, len3 = len(u5), len(u3)

			# Each strand's contiguous dot-run is a bulge/loop run
			if len5 > 0:
				bulge_runs.append(len5)
			if len3 > 0:
				bulge_runs.append(len3)

			# True bulge (unpaired on one strand only) vs internal loop (both strands)
			if (len5 > 0) ^ (len3 > 0):
				nr_true_bulge_bases += len5 + len3
			elif len5 > 0 and len3 > 0:
				nr_symmetric_loop_bases += len5 + len3

	if bulge_runs:
		max_len = max(bulge_runs)
		excess_runs = [r for r in bulge_runs if r > train_median_bulge_len]
		total_bulge_bases = sum(bulge_runs)
		excess_bulge_bases = sum(excess_runs)

		features["max_bulge_len"] = float(max_len)
		features["nr_excess_bulges"] = float(len(excess_runs))
		features["total_nr_bulge_bases"] = float(total_bulge_bases)
		features["pct_bulge"] = 100.0 * total_bulge_bases / n
		features["excess_bulge_nr_bases"] = float(excess_bulge_bases)
		features["pct_excess_bulges"] = 100.0 * excess_bulge_bases / n
		features["bulge_to_stem_ratio"] = (
			total_bulge_bases / nr_paired if nr_paired > 0 else 0.0
		)
		features["bulge_length_L2_norm"] = math.sqrt(sum(r * r for r in bulge_runs))
		features["nr_bases_true_bulges"] = float(nr_true_bulge_bases)
		features["nr_bases_symmetric_loops"] = float(nr_symmetric_loop_bases)
		features["pct_true_bulges"] = 100.0 * nr_true_bulge_bases / n
		features["pct_symmetric_loops"] = 100.0 * nr_symmetric_loop_bases / n

	return features


def _check_rod_features() -> None:
	"""Lock the ROD v2 definitions on a few hand-counted structures."""

	def expect(name: str, got: float, want: float) -> None:
		if abs(got - want) > 1e-9:
			raise ValueError("%s is %r, expected %r" % (name, got, want))

	def check(structure: str, sequence: str, expected: Dict[str, float]) -> None:
		got = extract_circular_rna_features(structure, sequence=sequence)
		expect("pS", got["pS"], got["pct_paired_bases"])
		total = got["pS"] + got["pH"] + got["pI"] + got["pB"] + got["pF"]
		if got["nr_paired_bases"] > 0:
			expect("state sum", total, 100.0)
		else:
			expect("state sum", total, 0.0)
		if got["nBranch"] == 0.0 and got["nr_paired_bases"] > 0:
			expect("pRod", got["pRod"], got["pS"])
			expect("pLongestRod", got["pLongestRod"], got["pS"])
		if got["nr_paired_bases"] > 0:
			if got["pLongestRod"] > got["pRod"] + 1e-9 or got["pRod"] > got["pS"] + 1e-9:
				raise ValueError("pLongestRod <= pRod <= pS failed for %s" % structure)
		for name, want in expected.items():
			expect(name, got[name], want)

	def same_circle(structures: List[str]) -> None:
		first = extract_circular_rna_features(structures[0])
		keys = (
			"pS",
			"pH",
			"pI",
			"pB",
			"pF",
			"pRod",
			"nBranch",
			"pLongestRod",
			"pOffRod",
			"pOnRod",
			"pSplitStem",
			"pLongestPairShare",
			"pJunction",
			"pDominantSpoke",
		)
		for structure in structures[1:]:
			got = extract_circular_rna_features(structure)
			for key in keys:
				expect(key, got[key], first[key])

	hairpin = "(((...)))"
	check(
		hairpin,
		"",
		{
			"pS": 100.0 * 6 / 9,
			"pH": 100.0 * 3 / 9,
			"pI": 0.0,
			"pB": 0.0,
			"pF": 0.0,
			"nBranch": 0.0,
			"pRod": 100.0 * 6 / 9,
			"pLongestRod": 100.0 * 6 / 9,
			"pOffRod": 0.0,
			"pOnRod": 100.0,
			"pSplitStem": 0.0,
			"pLongestPairShare": 100.0,
			"pJunction": 0.0,
			"pDominantSpoke": 100.0,
			"bKeto": 0.0,
			"uKeto": 0.0,
		},
	)
	check(
		hairpin,
		"GGGAAACCU",
		{"bKeto": 2.0 / 6.0, "uKeto": 1.0},
	)
	check(
		"(((.(...))))",
		"",
		{
			"pS": 100.0 * 8 / 12,
			"pH": 100.0 * 3 / 12,
			"pB": 100.0 * 1 / 12,
			"pI": 0.0,
			"pF": 0.0,
			"nBranch": 0.0,
		},
	)
	check(
		"((..((...))..))",
		"",
		{
			"pS": 100.0 * 8 / 15,
			"pH": 100.0 * 3 / 15,
			"pI": 100.0 * 4 / 15,
			"pB": 0.0,
			"pF": 0.0,
			"nBranch": 0.0,
		},
	)
	check(
		"(((...)))(((...)))(((...)))",
		"",
		{
			"pS": 100.0 * 18 / 27,
			"pH": 100.0 * 9 / 27,
			"pF": 0.0,
			"nBranch": 1.0,
			"pRod": 100.0 * 12 / 27,
			"pLongestRod": 100.0 * 6 / 27,
			"pOffRod": 100.0 * (18 - 12) / 18,
			"pOnRod": 100.0 * 12 / 18,
			"pSplitStem": 100.0 * (1.0 - 6.0 / 18.0),
			"pLongestPairShare": 100.0 * 6 / 18,
			"pJunction": 100.0 * (1.0 - 6.0 / 18.0),
			"pDominantSpoke": 100.0 * 6 / 18,
		},
	)
	check(
		"((...))......((...))",
		"",
		{
			"nBranch": 0.0,
			"pOffRod": 0.0,
			"pOnRod": 100.0,
			"pJunction": 0.0,
			"pDominantSpoke": 100.0,
			"pSplitStem": 50.0,
			"pLongestPairShare": 50.0,
		},
	)
	check(
		"((((...)))(((...))))",
		"",
		{
			"pS": 100.0 * 14 / 20,
			"pH": 100.0 * 6 / 20,
			"pF": 0.0,
			"nBranch": 1.0,
			"pRod": 100.0 * 12 / 20,
			"pLongestRod": 100.0 * 6 / 20,
		},
	)
	check(
		".....",
		"GGGGG",
		{
			"pS": 0.0,
			"pH": 0.0,
			"bKeto": 0.0,
			"uKeto": 1.0,
			"nBranch": 0.0,
			"pOffRod": 0.0,
			"pOnRod": 0.0,
			"pSplitStem": 0.0,
			"pLongestPairShare": 0.0,
			"pJunction": 0.0,
			"pDominantSpoke": 0.0,
		},
	)
	same_circle(
		[
			"...(((...)))",
			"(((...)))...",
			".(((...)))..",
			"..(((...))).",
		]
	)
	same_circle(
		[
			"(((...)))..(((...)))",
			"..(((...)))(((...)))",
			".(((...)))(((...))).",
		]
	)
	try:
		extract_circular_rna_features("(((...)))", sequence="GG")
	except ValueError:
		pass
	else:
		raise ValueError("a short sequence was accepted")


# --- Example Usage ---
if __name__ == "__main__":
	_check_rod_features()
	# # Structure with:
	# # - 2x2 internal loop: '..' (2 nt) on 5' strand and '..' (2 nt) on 3' strand
	# # - 1-nt true bulge: '.' (1 nt) on 5' strand only
	# sample_ss = "((..(((.(((...))))))..))"
	# feats = extract_circular_rna_features(
	# 	sample_ss, max_unpaired_gap=5, train_median_bulge_len=1.0
	# )

	# for k, v in feats.items():
	# 	print(f"{k:<42}: {v:.2f}")

	# zetavirus
	sample_ss = "........(((((..((((((((((((((((.((((((((((((((((.((...((....((.((((((((.(((((....))))).)))).)))).))....))...)).))))).))))))))))).)))).((((..(((.((((((((.(.((.(((((((.....))))))))).).))))))))..(((((((..........))))))).))))))).....((((((......)))))).)))))))).))))......)))))..........(((((((..((((.(((((((......)))))))...))))..)).))))).....((((((((...))))))))..(((((((((.((.(((((((((((.....)))))).))))).)))))))))))."
	feats = extract_circular_rna_features(
		sample_ss, max_unpaired_gap=5, train_median_bulge_len=1.0
	)
	for k, v in feats.items():
		print(f"{k:<42}: {v:.2f}")
