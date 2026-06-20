---
title: "Chapter 7 — Charge density & visualization"
sidebar_label: "7 · Charge density & viz"
sidebar_position: 8
---

# Chapter 7 — Charge density & visualization

Bands and DOS live in energy space. The **charge density** $n(\mathbf{r})$ lives
in *real* space — it is what DFT actually minimizes, and looking at it shows you
the bonds directly. This chapter uses `pp.x` to extract the GaAs valence density,
plots a slice through the Ga–As bonds, and points you at XCrySDen / VESTA for 3D.
It also closes out Phase 1 (ground state + electronic structure).

**Prerequisites:** [Chapter 6](./06-dos.md).

---

## 1. The quantity itself

The charge density is the sum of the occupied Kohn–Sham orbitals,

$$
n(\mathbf{r}) = \sum_{n}\int_\text{BZ}\frac{d\mathbf{k}}{V_\text{BZ}}\,
f_{n\mathbf{k}}\,|\psi_{n\mathbf{k}}(\mathbf{r})|^2 ,
$$

the self-consistent output of every SCF (it is literally what the loop in
Chapter 0 converged). It integrates to the number of valence electrons, and its
*shape* — where charge accumulates between atoms — is the picture of chemical
bonding.

## 2. `pp.x`: extract, then plot

Post-processing the density is a two-namelist job. `&inputpp` says **what** to
extract (via `plot_num`); `&plot` says **how** to slice and **what format** to
write:

```fortran
&inputpp
    prefix   = 'gaas'
    outdir   = './out'
    plot_num = 0            ! 0 = valence charge density
    filplot  = 'gaas.charge'
/
&plot
    iflag         = 2       ! 1=line, 2=plane, 3=3D
    output_format = 7       ! 7 = gnuplot 2D ; 5 = XSF (XCrySDen/VESTA) ; 6 = cube
    fileout       = 'gaas.rho_110.dat'
    e1(1)=1.0, e1(2)=1.0, e1(3)=0.0   ! first in-plane vector ([110])
    e2(1)=0.0, e2(2)=0.0, e2(3)=1.0   ! second in-plane vector ([001])
    x0(1)=0.0, x0(2)=0.0, x0(3)=0.0   ! origin of the plane (on a Ga)
    nx=180, ny=130
/
```

`plot_num` selects the field — a few of the most useful:

| `plot_num` | Field |
|---|---|
| 0 | valence charge density *n*(**r**) |
| 6 | spin density *n*↑ − *n*↓ (Chapter 9) |
| 7 | \|ψ\|² of selected band(s) / k-point(s) |
| 8 | electron localization function (ELF) |
| 9 | density **minus** superposition of atomic densities (deformation / bonding charge) |
| 11 | electrostatic potential (*V*<sub>bare</sub> + *V*<sub>H</sub>) |

`iflag` chooses the dimensionality (a 1D line, a 2D plane, or the full 3D grid),
and `output_format` chooses the file type — **`7`** for a 2D gnuplot file we can
read in Python, **`5`** for an `.xsf` you open in XCrySDen/VESTA, **`6`** for a
Gaussian `.cube`.

## 3. The workflow

```bash
cd code/07-charge
pw.x < gaas.scf.in     > scf.out         # self-consistent density
pp.x < gaas.pp_110.in  > pp_110.out      # 2D slice (1̄1̅0 plane) -> gaas.rho_110.dat
pp.x < gaas.pp_xsf.in  > pp_xsf.out      # full 3D grid       -> gaas.rho.xsf
```

We cut the **(1̄1̅0) plane** because it contains the Ga–As bonds: the bond runs
along [111], and this plane (spanned by [110] and [001]) holds whole
…Ga–As–Ga… chains.

## 4. Hands-on: the GaAs bond charge

[`notebooks/07-charge-density.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/07-charge-density.ipynb)
reads the 2D file and overlays the atom positions:

![GaAs valence charge density, (1̄1̅0) plane](/img/07-charge.png)

:::note The Ga 3d swamps the raw scale
Because this Ga pseudopotential keeps the **3d semicore** in valence, the raw
density spikes to ~1.5 e/bohr³ at the Ga sites — enough to wash out everything
else. We **cap the colour scale at 0.12** to reveal the chemistry. (Alternatively,
`plot_num = 9` plots the density *minus* the atomic superposition, isolating the
bonding charge directly — try it in the exercises.)
:::

With the scale capped, the physics is clear: charge does not just sit on the
atoms, it **bridges them** — a ridge of density along each Ga–As [111] bond, the
signature of a covalent bond. Look closely and the ridge is **pulled toward the
As**: As is more electronegative, so the shared charge is lopsided. That small
ionic asymmetry is exactly what makes GaAs *polar* — the origin of its Born
effective charges and the LO–TO phonon splitting we'll compute in Chapter 11.

## 5. Visualizing in 3D

The 2D slice is convenient for a figure, but real exploration is 3D. The
`gaas.rho.xsf` file (written with `output_format = 5`) opens directly in:

- **[XCrySDen](http://www.xcrysden.org/)** — `xcrysden --xsf gaas.rho.xsf`, then
  *Tools → Data Grid* to render an **isosurface**; drag the isovalue to watch the
  bonds appear.
- **[VESTA](https://jp-minerals.org/vesta/)** — *File → Open* the `.xsf`; VESTA
  shows the structure and lets you add a density isosurface or a 2D slice.

An isosurface at a small density value wraps the Ga–As bonds; raising it shrinks
the surface onto the atomic cores. The same `.xsf`/`.cube` route is how you'd
later visualize an ELF, a single Kohn–Sham orbital (`plot_num = 7`), or a spin
density.

## Exercises

1. **Count the electrons.** The density integrates to the number of valence
   electrons. How many should that be for GaAs, and which `pp.x` field would you
   integrate?
2. **Bonding charge directly.** Re-run with `plot_num = 9` (density minus atomic
   superposition) on the same plane. How does the picture change — where is
   charge *gained* and where is it *lost* relative to free atoms?
3. **A single orbital.** Use `plot_num = 7` with the highest valence band at Γ
   (`kpoint` and `kband`). Does its $|\psi|^2$ look bond-like or atom-like, and
   does that match the As-p character you found in the fatband/PDOS?
4. **ELF.** Compute `plot_num = 8` (ELF) on the bond plane. Where is the ELF
   close to 1, and what does that tell you about electron localization in the
   Ga–As bond?

<details>
<summary>Show solutions</summary>

**1. Count the electrons.** GaAs has **18** valence electrons (Ga 13 + As 5).
Integrating `plot_num = 0` (the valence charge density) over the cell returns 18
— it is the same number the SCF and the integrated DOS gave.

**2. Bonding charge (`plot_num = 9`).** Subtracting the free-atom densities
removes the big atomic peaks and shows the *redistribution*: charge is **gained**
in the bond region (the [111] ridge, shifted toward As) and **depleted** from the
back-bond/lobe regions around the atoms. This "deformation density" is the
cleanest way to see covalent bonding without the semicore swamping the scale.

**3. A single orbital (`plot_num = 7`).** The topmost valence state at Γ is
p-like and centred on As, with lobes pointing along the bonds — matching the
As-4p character from the PDOS and fatband. (A conduction state would instead look
s-like and more delocalized.)

**4. ELF.** The ELF approaches 1 in the **bond region between Ga and As**,
signalling a localized (paired) bonding electron — the hallmark of a covalent
bond. It is low in the interstitial regions. ELF is a popular complement to the
raw density precisely because it normalizes out the magnitude and highlights
*localization*.

</details>

---

## Phase 1 complete 🎉

You now have the full ground-state + electronic-structure toolkit: SCF,
plane-wave/pseudopotential and k-point convergence, relaxation, band structures,
DOS/PDOS, and real-space densities — all on one consistent GaAs example.

**Phase 2** moves to the physics that needs more than a plain PBE ground state:
[exchange–correlation functionals & DFT+U](./08-functionals.md) (Chapter 8), spin
and magnetism (Chapter 9), and lattice dynamics via DFPT — phonons and the LO–TO
splitting this chapter's bond charge hinted at (Chapters 10–11) — finishing with
running QE at scale on HPC (Chapter 12).

**Next:** [Chapter 8 — Exchange–correlation functionals & DFT+U](./08-functionals.md).
