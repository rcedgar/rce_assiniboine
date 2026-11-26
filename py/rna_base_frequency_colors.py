import numpy as np


def base_freqs_to_rgb(
	freqs: np.ndarray,
	alpha_gc: float = 2.2,
	alpha_skew: float = 3.5,
) -> np.ndarray:
	"""
	Maps nucleotide frequencies [f_A, f_U, f_G, f_C] (shape: (..., 4))
	to sRGB colors in [0, 1] with exact reverse-complement symmetry.

	- GC content controls perceptual lightness L (invariant under complement).
	- Purine/Pyrimidine skew (R-Y) controls OKLab 'a' (flips sign under complement).
	- Keto/Amino skew (K-M) controls OKLab 'b' (flips sign under complement).
	"""
	freqs = np.asarray(freqs, dtype=float)
	fA, fU, fG, fC = freqs[..., 0], freqs[..., 1], freqs[..., 2], freqs[..., 3]

	# 3 orthogonal simplex coordinates in [-1, 1]
	x_gc = (fG + fC) - (fA + fU)  # Invariant on reverse complement
	x_ry = (fA + fG) - (fU + fC)  # Negated on reverse complement
	x_km = (fG + fU) - (fA + fC)  # Negated on reverse complement

	# Non-linear stretch to separate biologically plausible compositions
	u_gc = np.tanh(alpha_gc * x_gc)
	u_ry = np.tanh(alpha_skew * x_ry)
	u_km = np.tanh(alpha_skew * x_km)

	# Map to OKLab coordinates (L in [0.42, 0.86], chroma radius <= 0.14)
	L = 0.64 + 0.22 * u_gc
	a = 0.13 * u_ry
	b = 0.13 * u_km

	# Convert OKLab -> linear sRGB
	l_ = L + 0.3963377774 * a + 0.2158037573 * b
	m_ = L - 0.1055613458 * a - 0.0638541728 * b
	s_ = L - 0.0894841775 * a - 1.2914855480 * b

	l = l_ ** 3
	m = m_ ** 3
	s = s_ ** 3

	r_lin = +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s
	g_lin = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s
	b_lin = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s

	rgb_lin = np.stack([r_lin, g_lin, b_lin], axis=-1)
	rgb_lin = np.clip(rgb_lin, 0.0, 1.0)

	# Gamma encode to standard sRGB
	srgb = np.where(
		rgb_lin <= 0.0031308,
		12.92 * rgb_lin,
		1.055 * (rgb_lin ** (1.0 / 2.4)) - 0.055,
	)
	return np.clip(srgb, 0.0, 1.0)


# --- Verification on Biologically Plausible Compositions ---
if __name__ == "__main__":
	examples = {
		"Uniform (25/25/25/25)":          [0.25, 0.25, 0.25, 0.25],
		"GC-rich (15/15/35/35)":          [0.15, 0.15, 0.35, 0.35],
		"AU-rich (35/35/15/15)":          [0.35, 0.35, 0.15, 0.15],
		"Purine-skewed (32/18/32/18)":    [0.32, 0.18, 0.32, 0.18],
		"Pyrimidine-skewed (18/32/18/32)":[0.18, 0.32, 0.18, 0.32],
		"G-biased loop (20/20/38/22)":    [0.20, 0.20, 0.38, 0.22],
		"C-biased comp (20/20/22/38)":    [0.20, 0.20, 0.22, 0.38],
	}

	for name, f in examples.items():
		rgb = base_freqs_to_rgb(np.array(f))
		hex_code = "#{:02x}{:02x}{:02x}".format(*(rgb * 255).astype(int))
		print(f"{name:<34} -> {hex_code}")