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
bulge_to_stem_ratio # total_nr_bulge_bases/
bulge_length_L2_norm # sqrt(sum of bulge length squared)
nr_bases_true_bulges # unpaired on on strand only
nr_bases_symmetric_loops # unpaired on both strands
pct_true_bulges # 100*nr_bases_true_bulges/nr_bases
pct_symmetric_loops # 100*nr_bases_symmetric_loops/nr_bases
'''

import math
from typing import Dict, List, Optional, Tuple


def extract_circular_rna_features(
	dot_bracket: str,
	max_unpaired_gap: int = 5,
	train_median_bulge_len: float = 1.0,
) -> Dict[str, float]:
	"""
	Parses a 3-symbol dot-bracket string ('(', ')', '.') representing a circular RNA
	and extracts global, strict-stem, composite-stem, and bulginess features.

	Args:
		dot_bracket: Circular secondary structure in 3-symbol dot-bracket notation.
		max_unpaired_gap: Max consecutive unpaired nucleotides allowed per strand
			when bridging helices into composite stems and detecting bulges/loops.
		train_median_bulge_len: Median bulge run length from the training set
			used to define excess bulges (> train_median_bulge_len).

	Returns:
		Dictionary containing all requested structural and bulginess features.
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


# --- Example Usage ---
if __name__ == "__main__":
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
