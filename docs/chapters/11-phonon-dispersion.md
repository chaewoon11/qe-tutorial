---
title: "Chapter 11 — Phonon dispersion"
sidebar_label: "11 · Phonon dispersion"
sidebar_position: 12
---

# Chapter 11 — Phonon dispersion

Chapter 10 gave the phonons at a single point (Γ). A **dispersion** $\omega_\nu(\mathbf{q})$
needs frequencies along a whole path through the Brillouin zone. We can't afford a
DFPT calculation at every q on the path, so we use the standard three-step trick:
DFPT on a coarse **q-grid** → real-space **force constants** (`q2r.x`) → cheap
**Fourier interpolation** to any q (`matdyn.x`). Along the way we impose the
acoustic sum rule and add the **LO–TO splitting** from the Born charges of
Chapter 10.

**Prerequisites:** [Chapter 10](./10-phonons-gamma.md) (DFPT, Born charges, ASR).

---

## 1. Why interpolate: q-grid → force constants → any q

The interatomic force constants $C(\mathbf{R})$ are **short-ranged** in real space
(forces between distant atoms die off). So the dynamical matrix $D(\mathbf{q})$ —
their Fourier transform — is *smooth* in q and fully determined by sampling it on
a modest, uniform **q-grid**. The procedure:

1. **`ph.x` with `ldisp=.true.`** computes $D(\mathbf{q})$ by DFPT on a uniform
   grid (here **4×4×4**), writing one dynamical matrix per irreducible q.
2. **`q2r.x`** inverse-Fourier-transforms those into real-space force constants
   $C(\mathbf{R})$ (file `gaas.fc`).
3. **`matdyn.x`** Fourier-transforms $C(\mathbf{R})$ *back* to **any** q you ask
   for — a dense path for the dispersion, or a dense grid for the phonon DOS — at
   negligible cost.

This is the phonon analogue of Wannier interpolation: do the expensive
first-principles work on a coarse grid, then interpolate freely.

## 2. Two physics fixes in the interpolation

**Acoustic sum rule (ASR).** As in Chapter 10, translational invariance requires
the three acoustic branches → 0 at Γ. Finite grids violate it slightly; `q2r.x`
(`zasr`) and `matdyn.x` (`asr`) re-impose it on the force constants so the
acoustic branches start cleanly from zero. (Recall: ASR only *enforces* this — it
relies on the underlying calculation already being well converged, which we
checked via the raw Γ acoustic frequency in Chapter 10.)

**LO–TO splitting.** In a **polar** crystal a *longitudinal* optical (LO) phonon
near Γ creates a macroscopic electric field that stiffens it above the
*transverse* optical (TO) mode — the **LO–TO splitting**. This is the
**non-analytic** part of the dynamical matrix at $\mathbf{q}\to0$, built from the
**Born effective charges** $Z^*$ and the dielectric constant $\varepsilon_\infty$
we computed with `epsil=.true.` in Chapter 10:

$$
D^{\text{NA}}_{I\alpha,J\beta}(\mathbf{q}\to0) \propto
\frac{(\mathbf{q}\cdot Z^*_I)_\alpha\,(\mathbf{q}\cdot Z^*_J)_\beta}
{\mathbf{q}\cdot\varepsilon_\infty\cdot\mathbf{q}} .
$$

`q2r.x` stores $Z^*$ and $\varepsilon_\infty$ in `gaas.fc`, and `matdyn.x` adds
this term automatically as the path approaches Γ — producing the LO frequency
above the TO.

## 3. The workflow

All of this runs from the converged SCF of Chapter 10 (`ecutwfc=70`,
`ecutrho=560`, 12×12×12 k, `conv_thr=1e-20`). The DFPT q-grid step is the
expensive one, so the whole chain ran on **Nurion's `debug` queue**:

```bash
pw.x    < gaas.scf.in        > gaas.scf.out        # converged density
ph.x    < gaas.ph.disp.in    > gaas.ph.disp.out    # DFPT on the 4x4x4 q-grid (the heavy step)
q2r.x   < gaas.q2r.in        > gaas.q2r.out        # -> gaas.fc (force constants + Z*, eps)
matdyn.x < gaas.matdyn.in    > gaas.matdyn.out     # -> gaas.freq.gp  (dispersion along L-G-X-W-K-G)
matdyn.x < gaas.matdyn.dos.in > gaas.matdyn.dos.out # -> gaas.phdos   (phonon DOS, 24^3 grid)
```

```fortran title="code/11-phonon-dispersion/gaas.ph.disp.in"
GaAs phonon dispersion: DFPT on a 4x4x4 q-grid
&inputph
    prefix   = 'gaas'
    outdir   = './out'
    fildyn   = 'gaas.dyn'
    ldisp    = .true.
    nq1 = 4, nq2 = 4, nq3 = 4
    tr2_ph   = 1.0d-16
    nmix_ph  = 12
    epsil    = .true.
    trans    = .true.
/
```

## 4. Hands-on: the GaAs phonon dispersion

[`notebooks/11-phonon-dispersion.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/11-phonon-dispersion.ipynb)
plots `gaas.freq.gp` along the path with the phonon DOS alongside:

![GaAs phonon dispersion and DOS](/img/11-phonon-dispersion.png)

Read it:

- **Acoustic branches** rise linearly from zero at Γ (their slopes are the sound
  velocities) and span the lower frequencies.
- **Optical branches** sit on top. At Γ they split: the **TO** at 269.7 cm⁻¹
  (doubly degenerate) and the **LO** at 284.0 cm⁻¹ — an **LO–TO splitting** of
  ~14 cm⁻¹, the direct consequence of the Born charges (a non-polar crystal like
  Si would have LO = TO at Γ).
- The **phonon DOS** peaks line up with the flat parts of the dispersion (van
  Hove singularities), exactly as the electronic DOS did with the bands.

The **TO at 269.7 cm⁻¹** matches experiment (~268) well. The **LO–TO splitting
(~14 cm⁻¹) is too small**, though — experiment is ~24 cm⁻¹. The splitting scales
as $(Z^*)^2/\varepsilon_\infty$, and PBE's over-large $\varepsilon_\infty$
(Chapter 10) suppresses it: the band-gap problem propagating one step further
into the lattice dynamics.

:::tip Convergence carries over
The dispersion is only as good as the q-grid (does the interpolation resolve the
real $C(\mathbf{R})$?) and the underlying DFPT convergence (Chapter 10). A 4×4×4
q-grid is fine for a 2-atom cell; denser cells / longer-ranged force constants
need denser q. And the LO–TO splitting is only as good as the $Z^*$ and
$\varepsilon_\infty$ — recall PBE over-polarizes GaAs (Chapter 10), so the
splitting here inherits that error.
:::

## Exercises

1. **Why a coarse q-grid suffices.** Explain, in terms of the *range* of the
   force constants $C(\mathbf{R})$, why DFPT on a 4×4×4 grid can reproduce the
   dispersion at thousands of q-points along the path. When would 4×4×4 *not* be
   enough?
2. **Spot the LO–TO splitting.** On the plot, the optical branches separate at Γ
   but the splitting *collapses* as you move away from Γ along some directions.
   Why is the non-analytic term a Γ-only (long-wavelength) effect?
3. **Non-polar contrast.** If you ran this whole workflow for silicon, what would
   be different at Γ for the optical modes, and why? (Think about $Z^*$.)
4. **DOS from dispersion.** The phonon DOS used a 24×24×24 q-grid via `matdyn.x`,
   not another DFPT run. Why is that cheap, and what sets the accuracy of the DOS?

<details>
<summary>Show solutions</summary>

**1. Why a coarse q-grid suffices.** $C(\mathbf{R})$ is short-ranged, so it is
nonzero only out to a few neighbour shells; a uniform q-grid of spacing
$2\pi/(N a)$ resolves force constants out to range $\sim N a$. Once the grid is
dense enough to capture the real range of $C(\mathbf{R})$, `matdyn.x` reproduces
$D(\mathbf{q})$ everywhere exactly. 4×4×4 fails when the force constants are
*long-ranged* — large/complex cells, or systems with long-range interactions
beyond the LO–TO term — which then need a denser q-grid.

**2. LO–TO is a Γ-only effect.** The non-analytic term depends on the *direction*
$\hat{\mathbf{q}}$ and survives only in the $\mathbf{q}\to0$ limit, where a
longitudinal optical displacement sets up a macroscopic, long-wavelength electric
field. At finite q that field is screened/curtailed, so the LO and TO branches
merge into the regular dispersion away from Γ.

**3. Non-polar contrast.** Silicon has two *identical* atoms by symmetry, so its
Born charges vanish ($Z^* = 0$): no macroscopic field, **no LO–TO splitting** —
the optical modes stay triply degenerate at Γ. The splitting is unique to polar
crystals like GaAs.

**4. DOS from dispersion.** Once `gaas.fc` exists, `matdyn.x` evaluates
frequencies on the 24³ grid by cheap Fourier interpolation — no new DFPT. The DOS
accuracy is set by (a) the interpolation grid density (24³ here) and (b) the
quality of `gaas.fc` itself, i.e. the original 4×4×4 DFPT q-grid and its
convergence.

</details>

---

**Next:** Chapter 12 — Running QE at scale *(coming soon)*: parallelization
levels (`-npool` / `-ndiag` / `-nband`), restarts, and HPC job scripts — the
practical layer behind the phonon jobs we just ran on Nurion.
