#!/usr/bin/env python3

import sys
import numpy as np


import numpy as np


def make_dna_grid(k: int = 8) -> tuple[np.ndarray, dict[tuple[int, int, int, int], int]]:
	"""
	Build a structured {A, C, G, T} grid with step 1/k and an O(1) lookup index.
	Choose k as a multiple of 4 (e.g., 4 -> 35 bins, 8 -> 165 bins, 12 -> 455 bins)
	so (0.25, 0.25, 0.25, 0.25) is an exact grid point.
	"""
	points = []
	index_map = {}
	for a in range(k + 1):
		for c in range(k + 1 - a):
			for g in range(k + 1 - a - c):
				t = k - a - c - g
				index_map[(a, c, g, t)] = len(points)
				points.append((a / k, c / k, g / k, t / k))
	return np.array(points, dtype=float), index_map


def lookup_bin(freqs: np.ndarray, k: int, index_map: dict) -> int:
	"""
	Snap any [p_A, p_C, p_G, p_T] vector to the nearest exact simplex grid bin in O(1).
	"""
	scaled = np.asarray(freqs, dtype=float) * k
	ints = np.rint(scaled).astype(int)
	diff = k - ints.sum()

	if diff != 0:
		# Adjust the element whose fractional part was closest to the rounding boundary
		residuals = scaled - ints
		order = np.argsort(residuals) if diff < 0 else np.argsort(-residuals)
		for idx in order[: abs(diff)]:
			ints[idx] += 1 if diff > 0 else -1

	return index_map[tuple(ints)]

def reject(g):
	assert len(g) == 4
	for x in g:
		if x > 0.99:
			return True
	return False


# Example with k=8 (165 bins—only slightly above 100, and hits 0.25 exactly)
k = 8
grid, idx_map = make_dna_grid(k)
N = len(grid)

L = 200
n = 0
print("# " + " ".join(sys.argv))
print("# L	fA	fC	fG	fT")
while L < 6000:
	L = int(L*1.2)
	grid, _ = make_dna_grid(k)
	assert len(grid) == N
	# print(grid.shape)   # (32, 4)
	# for i in range(N):
	for g in grid:
		if reject(g):
			continue
		s = str(L)
		for x in g:
			s += "\t%.3f" % x
		print(s)
		n += 1
print(f"# {k=} N={len(grid)}, {n=}")
