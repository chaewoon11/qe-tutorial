---
title: "Chapter 5 — Band structure"
sidebar_label: "5 · Band structure"
sidebar_position: 6
---

# Chapter 5 — Band structure

The SCF gave us occupied Kohn–Sham eigenvalues on a uniform grid. A **band
structure** plots those eigenvalues $\varepsilon_{n\mathbf{k}}$ along a path
through high-symmetry points of the Brillouin zone — the single most recognisable
picture in solid-state physics. This chapter builds the GaAs band structure,
reads off its **direct gap**, and meets the infamous DFT **band-gap problem**.

**Prerequisites:** [Chapter 4](./04-relaxation.md).

---

## 1. Why a separate, non-self-consistent step

The self-consistent density $n(\mathbf{r})$ only needs a uniform k-grid
(Chapter 3). But a band plot needs eigenvalues at *many* k-points along a path —
points that have tiny, unequal weights and would be a waste to include in the
SCF. The standard trick splits the work:

1. **`scf`** on a converged uniform grid → the self-consistent density and
   potential $v_\text{KS}[n]$.
2. A **non-self-consistent** step that *fixes* that potential and simply
   diagonalizes $\hat{H}_\text{KS}$ at each k-point along the path — no
   density update. For a band path this is `calculation = 'bands'` (a variant of
   `nscf` that expects an explicit k-path).

Because step 2 reuses the frozen potential, it is cheap, and you can throw as
many path k-points and empty bands (`nbnd`) at it as you like.

## 2. The k-path and `nbnd`

A zincblende (fcc) Brillouin zone has standard high-symmetry points
$\Gamma, X, L, W, K$. We walk the conventional path **L–Γ–X–W–K–Γ**, giving QE
the corner points in Cartesian $2\pi/a$ units (`tpiba_b`) and the number of
points to interpolate along each segment:

```fortran title="code/05-bands/gaas.bands.in (excerpt)"
&system
    nbnd = 16          ! 9 occupied + 7 empty: include conduction bands to see the gap
/
K_POINTS tpiba_b
6
  0.5  0.5  0.5  40    ! L
  0.0  0.0  0.0  40    ! Γ
  1.0  0.0  0.0  40    ! X
  1.0  0.5  0.0  20    ! W
  0.75 0.75 0.0  40    ! K
  0.0  0.0  0.0  1     ! Γ
```

The crucial change from the SCF input is **`nbnd = 16`**: the SCF only needed the
9 occupied bands, but to *see* a gap you must also compute empty (conduction)
bands.

## 3. The three-step workflow

```bash
cd code/05-bands
pw.x   < gaas.scf.in      > scf.out        # 1. self-consistent density
pw.x   < gaas.bands.in    > bands.out      # 2. eigenvalues along the path
bands.x < gaas.bands.pp.in > bands.pp.out  # 3. collect into a plottable file
```

The post-processor `bands.x` reads the raw eigenvalues and writes
`gaas.bands.dat.gnu` — a simple two-column (k-distance, energy) file, one band
per block — plus it prints the **x-coordinate of each high-symmetry point**,
which we use for the axis ticks:

```text
high-symmetry point:  0.5000 0.5000 0.5000   x coordinate   0.0000   (L)
high-symmetry point:  0.0000 0.0000 0.0000   x coordinate   0.8660   (Γ)
high-symmetry point:  1.0000 0.0000 0.0000   x coordinate   1.8660   (X)
...
```

## 4. Hands-on: the GaAs band structure

[`notebooks/05-band-structure.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/05-band-structure.ipynb)
parses the `.gnu` file, aligns the valence-band maximum to zero, and labels the
path:

![GaAs band structure (PBE)](/img/05-bands.png)

Read it from the bottom up:

- The flat five-band cluster near **−14.8 eV** is the **Ga 3d** semicore manifold
  — the same tightly-bound states we identified at Γ in Chapter 1, barely
  dispersing because they are core-like.
- The single band near **−12.7 eV** at Γ is the **As 4s** state.
- Three **p-like valence bands** rise to the **valence-band maximum at Γ**
  (degenerate there by symmetry).
- The lowest conduction band also bottoms out **at Γ** — so the smallest gap is
  **direct**.

Extracting the band edges (band 9 max, band 10 min):

```text
VBM = 7.618 eV (at Γ),  CBM = 8.132 eV (at Γ)
PBE band gap = 0.51 eV  (direct, at Γ)
```

GaAs is correctly predicted to be a **direct-gap** semiconductor — the physics
that makes it a laser and LED material.

## 5. The band-gap problem

That **0.51 eV** is badly wrong: the experimental room-temperature gap of GaAs is
**1.42 eV** ([Blakemore, *J. Appl. Phys.* **53**, R123
(1982)](https://doi.org/10.1063/1.331665)). The error is not a convergence
mistake — it is fundamental. **Kohn–Sham eigenvalue differences are not
quasiparticle (electron-removal/addition) energies**, and semilocal functionals
like PBE systematically *underestimate* band gaps, often by 50 % or more. GaAs is
a notorious case: with the Ga 3d in valence and PBE, the gap can come out
anywhere from a few tenths of an eV to nearly closed.

What the band structure *does* get right is the **shape**: the ordering of the
bands, their dispersions, the locations of the valley minima (Γ vs L vs X), and
the direct-vs-indirect character. Getting the gap *magnitude* right needs methods
beyond semilocal DFT — hybrid functionals (HSE), or many-body perturbation theory
(GW) — which are outside this tutorial's scope.

:::tip Trust shapes, distrust gap magnitudes
A PBE band structure is a reliable guide to band *topology* and *dispersions*,
but never quote its gap as a prediction. If you need an accurate gap, use a
hybrid or GW; if you only need effective masses or band ordering, PBE is often
fine.
:::

## Exercises

1. **Direct or indirect?** From the figure, confirm that both the VBM and the CBM
   lie at Γ. Where is the *second*-lowest conduction valley, and roughly how far
   above the Γ minimum is it? (GaAs's Γ–L–X valley ordering controls its
   transport.)
2. **Count the bands.** Why did we set `nbnd = 16`? What is the minimum `nbnd`
   that still shows the lowest conduction band, and what do you lose by using it?
3. **The gap error.** The PBE gap is 0.51 eV vs. 1.42 eV experiment. Express the
   error as a percentage. Why can't you fix it by tightening `ecutwfc`,
   `conv_thr`, or the k-path density?
4. **Path matters.** If you only computed bands along Γ–X (not through Γ
   properly), could you tell whether the gap is direct? What does that say about
   choosing a path that passes through *all* candidate valley positions?

<details>
<summary>Show solutions</summary>

**1. Direct or indirect.** Both edges are at Γ → direct gap. The next conduction
valley is at **L** (the conduction band dips again at the L end of the path),
roughly ~0.3–0.5 eV above the Γ minimum in PBE. This small Γ–L spacing is why
GaAs electrons can be scattered into the heavier L valley (the basis of the Gunn
effect).

**2. Count the bands.** 18 electrons → 9 occupied bands, so the SCF used 9. To
see the gap you need at least the 10th (first conduction) band, so the bare
minimum is `nbnd = 10`. We use 16 to show several conduction bands; using only 10
would plot just the lowest conduction band and miss the higher valleys and
band crossings.

**3. The gap error.** $(1.42-0.51)/1.42 \approx 64\,\%$ too small. It is a
*functional* error — the Kohn–Sham gap differs from the true quasiparticle gap by
the derivative discontinuity, which PBE lacks. No amount of basis, SCF, or
k-point convergence changes it; you need a different physical method (hybrid, GW).

**4. Path matters.** Along Γ–X alone you would see the Γ VBM and the Γ
conduction minimum, but you could not rule out a *lower* conduction minimum
elsewhere (e.g. at L). Determining direct vs. indirect requires sampling **all**
candidate extrema — which is exactly why the standard path visits Γ, X, L, K,
and W.

</details>

---

**Next:** Chapter 6 — Density of states & PDOS *(coming soon)*: from eigenvalues
to $g(E)$, and projecting onto atomic orbitals to see *what* the bands are made
of.
