A **perceptually uniform 3D color space (CIELAB or OKLCH)** mapped directly to the **3 canonical orthogonal basis vectors of a 4-base simplex** satisfies every constraint: exact 3-degree-of-freedom parameterization, exact reverse-complement symmetry, and compressed dynamic range so subtle differences in biologically plausible compositions (`30%–70%` GC, slight purine/pyrimidine skews) remain clearly distinguishable.

---

### 1. The Mathematical Parameterization (3 Orthogonal Coordinates)

Given nucleotide frequencies $(f_A, f_U, f_G, f_C)$ with $f_A + f_U + f_G + f_C = 1$ (swap $T$ for $U$ in DNA), decompose the 3D probability tetrahedron into three orthogonal axes:

1. **GC Content (Lightness / Luminance axis, $L$):**

$$x_{\text{GC}} = (f_G + f_C) - (f_A + f_U) \in [-1, 1]$$


* **Reverse-complement symmetry:** Under $A \leftrightarrow U$ and $G \leftrightarrow C$, $x_{\text{GC}}$ is **completely invariant** ($x_{\text{GC}}^{\text{rc}} = +x_{\text{GC}}$). Both strands share the exact same brightness/value.


2. **Purine vs. Pyrimidine Skew (Hue/Chroma axis 1, $a$):**

$$x_{\text{RY}} = (f_A + f_G) - (f_U + f_C) \in [-1, 1]$$


* **Reverse-complement symmetry:** Under complementation, purines $(A, G)$ swap with pyrimidines $(U, C)$, so $x_{\text{RY}}^{\text{rc}} = -x_{\text{RY}}$ (exact sign flip).


3. **Keto vs. Amino Skew (Hue/Chroma axis 2, $b$):**

$$x_{\text{KM}} = (f_G + f_U) - (f_A + f_C) \in [-1, 1]$$


* **Reverse-complement symmetry:** Under complementation, keto bases $(G, U)$ swap with amino bases $(A, C)$, so $x_{\text{KM}}^{\text{rc}} = -x_{\text{KM}}$ (exact sign flip).



---

### 2. Why This Mapping Works Geometrically

In **CIELAB** $(L^*, a^*, b^*)$ or **OKLab** $(L, a, b)$:

* **Uniform composition** $(0.25, 0.25, 0.25, 0.25)$ maps to $(x_{\text{GC}}, x_{\text{RY}}, x_{\text{KM}}) = (0, 0, 0)$, which lands on a **neutral mid-gray** ($L^* = 62, a^* = 0, b^* = 0$).
* **Complementary strands** share the exact same lightness $L^*$ and exact same saturation/chroma $C^* = \sqrt{a^{*2} + b^{*2}}$, while sitting **exactly $180^\circ$ opposite on the perceptual hue wheel** $(a^*, b^*) \leftrightarrow (-a^*, -b^*)$. Visually, complementary strands look like true perceptual opposites of equal visual weight.
* **Stretching the Plausible Biological Range:** Because extreme compositions (e.g., 95% A or 95% G) are rare, linear mapping wastes ~70% of the color gamut on empty corners. Passing each coordinate through a **hyperbolic tangent (`tanh`)** or **power-law stretch** concentrates color resolution inside the biologically realistic window ($f_i \in [0.15, 0.40]$):

$$u_{\text{GC}} = \tanh(\alpha_{\text{GC}} \cdot x_{\text{GC}}), \quad u_{\text{RY}} = \tanh(\alpha_{\text{skew}} \cdot x_{\text{RY}}), \quad u_{\text{KM}} = \tanh(\alpha_{\text{skew}} \cdot x_{\text{KM}})$$

With $\alpha_{\text{GC}} \approx 2.2$ and $\alpha_{\text{skew}} \approx 3.5$:

* A modest shift from **40% GC to 60% GC** ($x_{\text{GC}}: -0.2 \to +0.2$) spans nearly half the lightness range.
* A subtle **5% strand skew** (e.g., $A=0.30, U=0.20, G=0.25, C=0.25 \implies x_{\text{RY}} = +0.10$) already produces a distinct, readable tint without saturating out at extreme skews.

---

### 3. Perceptual Color Wheel Assignment

Mapping $(u_{\text{GC}}, u_{\text{RY}}, u_{\text{KM}})$ into **OKLab** (or CIELAB) gives four intuitive corner tints around neutral gray:

| Composition Bias | Dominant Bases | Coordinates $(x_{\text{GC}}, x_{\text{RY}}, x_{\text{KM}})$ | Resulting Visual Appearance |
| --- | --- | --- | --- |
| **Balanced / Uniform** | $A=U=G=C=0.25$ | $(0, 0, 0)$ | **Neutral Slate Gray** |
| **GC-Rich vs. AU-Rich** | High $G+C$ vs. High $A+U$ | $(+, 0, 0)$ vs. $(-, 0, 0)$ | **Bright / Light** vs. **Deep / Dark** |
| **Purine-Rich (R)** | $A + G$ enriched | $(0, +, 0)$ | **Warm Coral / Terracotta** ($+a$) |
| **Pyrimidine-Rich (Y)** | $U + C$ enriched | $(0, -, 0)$ | **Cool Teal / Cyan** ($-a$) |
| **Keto-Rich (K)** | $G + U$ enriched | $(0, 0, +)$ | **Ochre / Amber Gold** ($+b$) |
| **Amino-Rich (M)** | $A + C$ enriched | $(0, 0, -)$ | **Periwinkle / Violet-Blue** ($-b$) |

Because **Purine vs. Pyrimidine** ($a$-axis) and **Keto vs. Amino** ($b$-axis) are orthogonal:

* **A-rich** ($+R, -K$) $\rightarrow$ **Rose / Magenta**
* **U-rich** ($-R, +K$) $\rightarrow$ **Sage / Olive Green** (exact complement to A-rich)
* **G-rich** ($+R, +K$) $\rightarrow$ **Warm Orange / Amber**
* **C-rich** ($-R, -K$) $\rightarrow$ **Cerulean / Cobalt Blue** (exact complement to G-rich)

---

### 4. Ready-to-Use Python Implementation (OKLab $\rightarrow$ sRGB)

Using **OKLab** (Björn Ottosson's perceptual color space) avoids the hue-twisting artifacts of CIELAB in the blue/purple region and clamps cleanly into sRGB:

```python
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

```

### Tuning Tips for Your Dataset

* **Adjusting sensitivity:** Set `alpha_gc` and `alpha_skew` to $\approx \frac{1}{\text{IQR}(x)}$ of each coordinate in your dataset so that the middle 50% of your sequences span `[-0.46, +0.46]` in `tanh` space.
* **Colorblind accessibility:** If red–green discrimination (`a`-axis) is a concern for figures, rotate the $(u_{\text{RY}}, u_{\text{KM}})$ plane by $45^\circ$ before assigning $(a, b)$ so that the two skew axes align with **blue–yellow** and **magenta–cyan** instead of red–green—preserving exact $180^\circ$ complement symmetry while staying safe for deuteranopia/protanopia.