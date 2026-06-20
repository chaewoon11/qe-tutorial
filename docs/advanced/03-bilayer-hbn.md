---
id: bilayer-hbn
title: "A3 — Bilayer hBN: vdW, bands & 2D phonons"
sidebar_label: "A3 · Bilayer hBN (vdW, LO–TO)"
sidebar_position: 3
---

# A3 — Bilayer hBN: vdW, bands & 2D phonons

Graphene (A1–A2) is a single, non-polar sheet. **Bilayer hexagonal boron
nitride** adds the two ingredients graphene can't show: **van der Waals
interlayer bonding** (you need a dispersion correction) and a **polar 2D lattice**
(so the phonons have a genuinely 2D **LO–TO splitting**). It is also a **wide-gap
insulator** — the opposite electronic extreme from graphene's Dirac cone.

:::note Reference data
The DFPT phonon run here (a full 12×12×1 q-grid at `ecutwfc = 98 Ry` on a 4-atom
cell) is an **HPC-grade** calculation; the band and phonon data shown are
pre-computed reference results. The inputs below are complete and runnable.
:::

**Prerequisites:** [A1](./01-2d-cell.md), [A2](./02-dirac-bands.md), and the core
[bands](../chapters/bands.md) + [phonon-dispersion](../chapters/phonon-dispersion.md)
chapters.

---

## 1. Two layers held by van der Waals

hBN is isostructural to graphene — a honeycomb of alternating **B** and **N** —
but the two sublattices are *different atoms*, so it is polar and insulating.
A bilayer is two such sheets; the common **AA′** stacking puts B over N and N
over B. The layers are bound only by **van der Waals** forces, which semilocal
PBE *does not describe* — so we add an explicit dispersion correction:

```fortran
&system
    vdw_corr        = 'dft-d3'      ! Grimme D3 dispersion — binds the layers
    dftd3_threebody = .false.       ! 3-body term is incompatible with DFPT; turn it off
/
```

Without it the interlayer spacing is badly overestimated (PBE has no long-range
attraction). The `dftd3_threebody = .false.` is essential here because we follow
up with DFPT phonons, which the 3-body D3 term does not support.

The complete self-consistent input (the AA′ relaxed bilayer, 2D-isolated):

```fortran title="code/advanced/03-bilayer-hbn/bn.scf.in"
&control
    calculation     = 'scf'
    prefix          = 'bn'
    pseudo_dir      = '../../pseudos'
    outdir          = './out'
    verbosity       = 'low'
    wf_collect      = .true.
/
&system
    ibrav           = 0
    nat             = 4
    ntyp            = 2
    ecutwfc         = 98
    vdw_corr        = 'dft-d3'
    dftd3_threebody = .false.
    assume_isolated = '2D'
    occupations     = 'fixed'
/
&electrons
    diagonalization = 'david'
    mixing_beta     = 0.7
    conv_thr        = 1.0d-12
/
ATOMIC_SPECIES
  B  10.811   B.pbe-n-kjpaw_psl.1.0.0.UPF
  N  14.0067  N.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
   2.508918318   0.000000000   0.000000000
  -1.254459159   2.172786999   0.000000000
   0.000000000   0.000000000  45.000000000
ATOMIC_POSITIONS crystal
  B   0.3333333333   0.6666666667  -0.0380964600
  N   0.6666666667   0.3333333333  -0.0381138600
  N   0.3333333333   0.6666666667   0.0381138600
  B   0.6666666667   0.3333333333   0.0380964600
K_POINTS automatic
  12 12 1 0 0 0
```

## 2. Bands: a wide-gap insulator

Same `scf → bands → bands.x` recipe as A2, along Γ–M–K–Γ:

```fortran title="code/advanced/03-bilayer-hbn/bn.bands.in"
&control
    calculation     = 'bands'
    prefix          = 'bn'
    pseudo_dir      = '../../pseudos'
    outdir          = './out'
    verbosity       = 'low'
/
&system
    ibrav           = 0
    nat             = 4
    ntyp            = 2
    nbnd            = 10
    ecutwfc         = 98
    vdw_corr        = 'dft-d3'
    dftd3_threebody = .false.
    assume_isolated = '2D'
    occupations     = 'fixed'
/
&electrons
    diago_full_acc  = .true.
    conv_thr        = 1.0d-16
/
ATOMIC_SPECIES
  B  10.811   B.pbe-n-kjpaw_psl.1.0.0.UPF
  N  14.0067  N.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
   2.508918318   0.000000000   0.000000000
  -1.254459159   2.172786999   0.000000000
   0.000000000   0.000000000  45.000000000
ATOMIC_POSITIONS crystal
  B   0.3333333333   0.6666666667  -0.0380964600
  N   0.6666666667   0.3333333333  -0.0381138600
  N   0.3333333333   0.6666666667   0.0381138600
  B   0.6666666667   0.3333333333   0.0380964600
K_POINTS crystal_b
4
  0.0000000000   0.0000000000   0.0   80   ! Γ
  0.5000000000   0.0000000000   0.0   50   ! M
  0.3333333333   0.3333333333   0.0  100   ! K
  0.0000000000   0.0000000000   0.0    1   ! Γ
```

```fortran title="code/advanced/03-bilayer-hbn/bn.bands.pp.in"
&bands
    prefix  = 'bn'
    outdir  = './out'
    filband = 'bn.bands'
/
```

## 3. 2D phonons with `loto_2d`

The DFPT chain is exactly core [Chapter 11](../chapters/phonon-dispersion.md) —
`ph.x` (with `epsil` for the Born charges) → `q2r.x` → `matdyn.x` — but with one
2D-specific switch: **`loto_2d = .true.`** in *both* `q2r.x` and `matdyn.x`. It
replaces the 3D dipole–dipole long-range term with its **2D** form, giving the
correct long-wavelength LO–TO behaviour for a single (or few-layer) sheet.

```fortran title="code/advanced/03-bilayer-hbn/bn.ph.in"
&inputph
  prefix   = 'bn'
  outdir   = './out'
  fildyn   = 'bn.dyn'
  fildvscf = 'dvscf'
  ldisp    = .true.
  epsil    = .true.            ! polar -> Born charges (needed for LO-TO)
  nmix_ph  = 12
  tr2_ph   = 1.0d-18
  nq1 = 12, nq2 = 12, nq3 = 1
/
```

```fortran title="code/advanced/03-bilayer-hbn/bn.q2r.in"
&input
  fildyn  = 'bn.dyn'
  flfrc   = 'bn.fc'
  zasr    = 'simple'
  loto_2d = .true.            ! 2D long-range dipole term
/
```

```fortran title="code/advanced/03-bilayer-hbn/bn.matdyn.in"
&input
  flfrc   = 'bn.fc'
  asr     = 'simple'
  flfrq   = 'bn.freq'
  loto_2d = .true.
  q_in_band_form = .true.
  q_in_cryst_coord = .false.
/
4
0.0000000000  0.0000000000  0.0   80
0.0000000000  0.5773502692  0.0   50
0.3333333333  0.5773502692  0.0  100
0.0000000000  0.0000000000  0.0    1
```

## 4. Hands-on: bands and phonons together

[`notebooks/A3-bilayer-hbn.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/A3-bilayer-hbn.ipynb)
stacks both over the same Γ–M–K–Γ path:

![Bilayer hBN bands and phonons](/img/A3-bilayer-hbn.png)

- **Electronic structure (top).** Bilayer hBN is a **wide-gap insulator** — PBE
  gives an **indirect** gap of **4.53 eV** (VBM near K, CBM near M). That is the
  opposite extreme from graphene: same honeycomb, but B≠N opens a large gap.
- **Phonons (bottom).** The acoustic branches rise from Γ (the lowest is the
  out-of-plane **ZA flexural** mode). The high optical branches near
  **1350–1400 cm⁻¹** are the in-plane **E$_{2g}$** modes (hBN's Raman signature).
- **2D LO–TO.** Because hBN is polar and we used `loto_2d`, the LO branch shows
  the **2D** form of LO–TO splitting: it is **degenerate with TO exactly at Γ**
  but climbs **linearly** in $q$ as you leave Γ (here from ~1352 to ~1410 cm⁻¹
  just off Γ). In 3D the LO–TO split is a finite *jump* at Γ; in 2D the long-range
  dipole field is weaker and the split is **linear-in-q** — which is exactly what
  `loto_2d` captures and a plain 3D treatment would get wrong.

This is the payoff of using a **polar** 2D crystal: every 2D-specific feature —
vdW binding, a real gap, and the 2D LO–TO — shows up, none of which graphene
could provide.

## Exercises

1. **Why `dft-d3`?** Re-relax the bilayer *without* `vdw_corr`. What happens to
   the interlayer distance, and why does PBE alone fail for layered materials?
2. **2D vs 3D LO–TO.** Re-run `matdyn.x` with `loto_2d = .false.`. How does the
   behaviour of the LO branch near Γ change, and which is correct for a bilayer?
3. **Graphene vs hBN.** Both are honeycombs with the same path. Explain — in terms
   of the two sublattice atoms — why graphene is a zero-gap semimetal while hBN is
   a wide-gap insulator with LO–TO splitting.
4. **`threebody`.** Why must `dftd3_threebody = .false.` be set when the run is
   followed by DFPT? What would break otherwise?

<details>
<summary>Show solutions</summary>

**1. Why `dft-d3`.** Without a dispersion correction PBE has no long-range
attraction, so the layers drift apart — the interlayer spacing comes out far too
large (or unbound). vdW forces come from correlated charge fluctuations that
semilocal functionals miss; D3 adds them back as a pairwise $-C_6/r^6$ term.

**2. 2D vs 3D LO–TO.** With `loto_2d=.false.` the code applies the 3D
dipole–dipole term, giving a spurious *finite jump* of the LO mode at Γ. For a
genuine 2D sheet that is wrong: the correct behaviour is the **linear-in-q** rise
that `loto_2d=.true.` produces.

**3. Graphene vs hBN.** In graphene the two honeycomb sites are identical C atoms,
so by symmetry the π and π\* bands are degenerate at K (Dirac point) → zero gap,
and the Born charges vanish (no LO–TO). In hBN the sites are B and N with very
different electronegativity: this breaks the sublattice symmetry, opens a large
gap, and gives nonzero Born charges → LO–TO splitting.

**4. `threebody`.** DFPT (phonons) needs analytic second derivatives of the
energy; the 3-body Axilrod–Teller–Muto term in D3 is not implemented in the DFPT
code path, so `ph.x` would stop (or be inconsistent). Setting
`dftd3_threebody=.false.` keeps the (pairwise) dispersion that *is* supported.

</details>

---

That closes the **2D-materials track**: graphene's setup and Dirac bands
(A1–A2), and bilayer hBN's vdW binding, wide gap, and 2D LO–TO phonons (A3) —
the 2D-specific pieces layered on top of the core QE workflow.
