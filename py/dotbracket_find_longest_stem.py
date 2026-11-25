from typing import Dict, List, Optional, Tuple


def parse_circular_longest_stem(
	dot_bracket: str,
	max_unpaired_gap: int = 2,
	include_hairpin_loop: bool = False,
) -> Dict:
	"""
	Parses a 3-symbol dot-bracket string for a circular RNA, identifies the
	longest composite stem (bridging bulges/internal loops up to `max_unpaired_gap`),
	and annotates both paired and unpaired residues belonging to that stem.
	"""
	n = len(dot_bracket)
	if n == 0:
		return {"longest_stem_length": 0, "longest_stem_pairs": [], "annotations": []}

	# 1. Parse base pairs using a stack
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

	canonical_pairs: List[Tuple[int, int]] = sorted(
		(i, j) for i, j in pair_map.items() if i < j
	)

	if not canonical_pairs:
		return {
			"longest_stem_length": 0,
			"longest_stem_pairs": [],
			"annotations": [
				{
					"index": i,
					"symbol": c,
					"is_paired": False,
					"partner": None,
					"role": "unpaired_exterior",
					"in_longest_stem": False,
				}
				for i, c in enumerate(dot_bracket)
			],
		}

	# Helper: collect circular indices strictly between start and end moving +1
	def circular_between(start: int, end: int) -> List[int]:
		indices = []
		curr = (start + 1) % n
		while curr != end:
			indices.append(curr)
			curr = (curr + 1) % n
		return indices

	# Helper: find next inward pair and any unpaired bulge/internal-loop indices bridged
	def next_inward_step(
		i: int, j: int
	) -> Optional[Tuple[Tuple[int, int], List[int]]]:
		for step_5 in range(1, max_unpaired_gap + 2):
			ni = (i + step_5) % n
			if dot_bracket[ni] == ".":
				continue
			for step_3 in range(1, max_unpaired_gap + 2):
				nj = (j - step_3) % n
				if dot_bracket[nj] == ".":
					continue
				if pair_map.get(ni) == nj and ni != j and nj != i:
					unpaired_5 = circular_between(i, ni)
					unpaired_3 = circular_between(nj, j)
					return (ni, nj), unpaired_5 + unpaired_3
				return None
			return None
		return None

	# Helper: find previous outward pair
	def prev_outward_pair(i: int, j: int) -> Optional[Tuple[int, int]]:
		for step_5 in range(1, max_unpaired_gap + 2):
			pi = (i - step_5) % n
			if dot_bracket[pi] == ".":
				continue
			for step_3 in range(1, max_unpaired_gap + 2):
				pj = (j + step_3) % n
				if dot_bracket[pj] == ".":
					continue
				if pair_map.get(pi) == pj and pi != j and pj != i:
					return (pi, pj)
				return None
			return None
		return None

	# 2. Group pairs and bridged unpaired residues into composite stems
	visited = set()
	stems: List[Dict] = []

	for i, j in canonical_pairs:
		if (i, j) in visited or (j, i) in visited:
			continue

		# Walk outward to the outermost base pair of this stem
		curr = (i, j)
		while True:
			prev_pair = prev_outward_pair(*curr)
			if prev_pair and prev_pair not in visited and prev_pair != (i, j):
				curr = prev_pair
			else:
				break

		# Walk inward to collect pairs and internal unpaired residues
		stem_pairs: List[Tuple[int, int]] = []
		stem_internal_unpaired: List[int] = []

		while curr and curr not in visited:
			visited.add(curr)
			visited.add((curr[1], curr[0]))
			stem_pairs.append(curr)

			step_result = next_inward_step(*curr)
			if step_result:
				next_pair, bridged_unpaired = step_result
				stem_internal_unpaired.extend(bridged_unpaired)
				curr = next_pair
			else:
				curr = None

		# Check if the innermost pair caps an apical hairpin loop (all dots inside)
		inner_i, inner_j = stem_pairs[-1]
		enclosed = circular_between(inner_i, inner_j)
		hairpin_unpaired = (
			enclosed if (enclosed and all(dot_bracket[k] == "." for k in enclosed)) else []
		)

		stems.append({
			"pairs": stem_pairs,
			"internal_unpaired": stem_internal_unpaired,
			"hairpin_unpaired": hairpin_unpaired,
		})

	# 3. Identify the longest stem (ranked by number of base pairs)
	longest = max(stems, key=lambda s: len(s["pairs"]))

	paired_in_longest = {idx for pair in longest["pairs"] for idx in pair}
	internal_in_longest = set(longest["internal_unpaired"])
	hairpin_in_longest = set(longest["hairpin_unpaired"])

	# Classify all unpaired bases across the structure for clear role labels
	all_internal = {idx for s in stems for idx in s["internal_unpaired"]}
	all_hairpin = {idx for s in stems for idx in s["hairpin_unpaired"]}

	# 4. Build per-residue annotation
	annotations = []
	for idx, char in enumerate(dot_bracket):
		partner = pair_map.get(idx, None)
		is_paired = partner is not None

		if is_paired:
			role = "paired_stem"
			in_longest = idx in paired_in_longest
		elif idx in all_internal:
			role = "unpaired_bulge_or_internal_loop"
			in_longest = idx in internal_in_longest
		elif idx in all_hairpin:
			role = "unpaired_hairpin_loop"
			in_longest = include_hairpin_loop and (idx in hairpin_in_longest)
		else:
			role = "unpaired_junction_or_exterior"
			in_longest = False

		annotations.append({
			"index": idx,
			"symbol": char,
			"is_paired": is_paired,
			"partner": partner,
			"role": role,
			"in_longest_stem": in_longest,
		})

	return {
		"longest_stem_length": len(longest["pairs"]),
		"longest_stem_pairs": longest["pairs"],
		"longest_stem_unpaired": sorted(longest["internal_unpaired"]),
		"annotations": annotations,
	}


# --- Example Usage ---
circular_db = "((..(((...)))..))"
result = parse_circular_longest_stem(
	circular_db, max_unpaired_gap=2, include_hairpin_loop=False
)

print(f"Structure:              {circular_db}")
print(f"Longest stem (bp):      {result['longest_stem_length']}")
print(f"Stem base pairs:        {result['longest_stem_pairs']}")
print(f"Stem unpaired residues: {result['longest_stem_unpaired']}\n")

print(f"{'Idx':<4} {'Sym':<4} {'Paired':<7} {'Partner':<8} {'Role':<32} {'In Longest Stem'}")
print("-" * 76)
for res in result["annotations"]:
	partner_str = str(res["partner"]) if res["partner"] is not None else "-"
	stem_mark = "YES" if res["in_longest_stem"] else ""
	print(
		f"{res['index']:<4} {res['symbol']:<4} {str(res['is_paired']):<7} "
		f"{partner_str:<8} {res['role']:<32} {stem_mark}"
	)