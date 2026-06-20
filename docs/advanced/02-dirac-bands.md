---
id: dirac-bands
title: "A2 — Dirac bands of graphene"
sidebar_label: "A2 · Dirac bands"
sidebar_position: 2
---

# A2 — Dirac bands of graphene

With the converged, relaxed graphene cell from [A1](./01-2d-cell.md), we compute
its band structure. Graphene's bands are the reason it is famous: the valence and
conduction π bands meet at a single point — the **Dirac cone** — making it a
zero-gap **semimetal** with massless-like carriers.

**Prerequisites:** [A1](./01-2d-cell.md) (relaxed cell) and the core
[band-structure chapter](../chapters/bands.md) (scf → bands → `bands.x`).

---

## 1. What to expect: σ, π, and the Dirac point

Each carbon contributes four valence electrons. Three go into in-plane **sp²
σ bonds** (deep, strongly-bound bands); the fourth goes into the out-of-plane
**π** orbital. A tight-binding model of the π electrons on the honeycomb lattice
gives two bands that **touch at the corner K** of the hexagonal Brillouin zone,
and because there is exactly one π electron per atom the **Fermi level sits right
at that touching point**. Near K the bands are **linear** in $\mathbf{k}$,

$$
E(\mathbf{k}) \approx \pm \hbar v_F\,\lvert \mathbf{k} - \mathbf{K}\rvert ,
$$

a **Dirac cone** — low-energy electrons behave like massless particles moving at
the Fermi velocity $v_F$. We'll see all of this fall straight out of DFT.

## 2. Workflow

Same three steps as the core band chapter — `scf` → `bands` → `bands.x` — on the
**A1 relaxed cell**, with the 2D-track thresholds (`scf` 1e-20, `bands` 1e-16):

```bash
cd code/advanced/02-dirac-bands
pw.x    < graphene.scf.in      > scf.out        # density (note the Fermi energy)
pw.x    < graphene.bands.in    > bands.out      # eigenvalues along Γ–M–K–Γ
bands.x < graphene.bands.pp.in > bands.pp.out   # -> graphene.bands.dat.gnu
```

The self-consistent step is the converged, relaxed monolayer:

```fortran title="code/advanced/02-dirac-bands/graphene.scf.in"
&control
    calculation = 'scf'
    prefix      = 'graphene'
    outdir      = './out'
    pseudo_dir  = '../../pseudos'
    verbosity   = 'high'
/
&system
    ibrav = 0
    nat = 2
    ntyp = 1
    ecutwfc = 50.0
    ecutrho = 400.0
    assume_isolated = '2D'
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.01
/
&electrons
    conv_thr = 1.0d-20
    mixing_beta = 0.7
/
ATOMIC_SPECIES
  C 12.011 C.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
   2.467137   0.000000   0.000000
  -1.233568   2.136488   0.000000
   0.000000   0.000000  20.000000
ATOMIC_POSITIONS crystal
  C 0.000000000 0.000000000 0.0
  C 0.333333333 0.666666667 0.0
K_POINTS automatic
  12 12 1 0 0 0
```

The band step uses `calculation='bands'` with the hexagonal path in
**`crystal_b`** coordinates (Γ, M $=(\tfrac12,0)$, K $=(\tfrac13,\tfrac13)$) and
enough bands to show the conduction states:

```fortran title="code/advanced/02-dirac-bands/graphene.bands.in"
&control
    calculation = 'bands'
    prefix      = 'graphene'
    outdir      = './out'
    pseudo_dir  = '../../pseudos'
/
&system
    ibrav = 0
    nat = 2
    ntyp = 1
    ecutwfc = 50.0
    ecutrho = 400.0
    assume_isolated = '2D'
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.01
    nbnd = 12
/
&electrons
    conv_thr = 1.0d-16
/
ATOMIC_SPECIES
  C 12.011 C.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
   2.467137   0.000000   0.000000
  -1.233568   2.136488   0.000000
   0.000000   0.000000  20.000000
ATOMIC_POSITIONS crystal
  C 0.000000000 0.000000000 0.0
  C 0.333333333 0.666666667 0.0
K_POINTS crystal_b
4
  0.0000 0.0000 0.0  60    ! Γ
  0.5000 0.0000 0.0  40    ! M
  0.3333 0.3333 0.0  60    ! K
  0.0000 0.0000 0.0   1    ! Γ
```

```fortran title="code/advanced/02-dirac-bands/graphene.bands.pp.in"
&bands
    prefix = 'graphene'
    outdir = './out'
    filband = 'graphene.bands.dat'
/
```

## 3. Hands-on: the Dirac cone

The SCF reports a Fermi level of **−4.225 eV**; we align the bands to it.
[`notebooks/A2-dirac-bands.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/A2-dirac-bands.ipynb)
plots the full path and a zoom on K:

![Graphene Dirac bands](/img/A2-dirac-bands.png)

The physics is exactly the tight-binding prediction:

- The deep parabolic bands (−8 to −20 eV) are the **σ** (sp²) framework.
- The **π / π\*** bands cross at **K**, precisely at the Fermi level — the
  **Dirac point** (our computed gap there is ~0.002 eV, i.e. zero to numerical
  precision).
- Zoomed in, the bands are **linear** around K: the Dirac cone.

Graphene is correctly a **semimetal** — no gap, but the density of states
vanishes *at* the Dirac point. (No band-gap problem here: there is no gap to
underestimate. PBE nails graphene's topology and Fermi velocity reasonably well.)

## Exercises

1. **Count the bands.** Graphene has 8 valence electrons per cell. How many
   occupied bands is that, and which one is the π band at K? Why does the Fermi
   level land exactly at the Dirac point (not above or below)?
2. **Miss the cone.** Re-run `bands.in` with a path that goes Γ–M–Γ only (no K).
   Could you tell graphene is a semimetal from that? What does this say about
   choosing a k-path?
3. **Fermi velocity.** From the zoom, estimate the slope $dE/dk$ of the Dirac
   cone (convert the x-axis from 2π/a units to Å⁻¹). How does your $v_F$ compare
   to the textbook ~10⁶ m/s?
4. **Why `crystal_b`?** We gave K as $(\tfrac13,\tfrac13,0)$ in `crystal_b`. Why
   are crystal coordinates more convenient than Cartesian (`tpiba_b`) for a
   hexagonal cell?

<details>
<summary>Show solutions</summary>

**1. Count the bands.** 8 electrons → **4 occupied bands** (spin-degenerate):
three σ and one π. At K the π band touches the empty π\*. The Fermi level sits at
the Dirac point because the π band holds exactly one electron per atom (half-
filling of the two π states) — neutral graphene fills states up to exactly the
touching point.

**2. Miss the cone.** Along Γ–M–Γ you never visit K, so you'd see an apparent gap
and wrongly conclude graphene is an insulator/semiconductor. The Dirac physics
lives *only* at K — the k-path **must** pass through it. Always route the path
through every symmetry point where the interesting physics could be.

**3. Fermi velocity.** Read $\Delta E/\Delta k$ from the linear region near K
(convert $k$ to Å⁻¹ using $2\pi/a$ with $a=2.467$ Å). The slope gives
$\hbar v_F$, and $v_F \approx 0.8\text{–}1\times10^6$ m/s — close to the
experimental value (PBE slightly underestimates it).

**4. Why `crystal_b`.** In crystal coordinates the high-symmetry points have
simple rational values fixed by symmetry (K is *always* $(\tfrac13,\tfrac13)$),
independent of the lattice constant. In Cartesian `tpiba_b` you'd have to work out
$(\tfrac13,\tfrac13)$ in $2\pi/a$ units (involving $\sqrt3$ factors) — easy to get
wrong for a hexagonal cell.

</details>

---

**Next:** A3 — 2D phonons of graphene *(coming soon)*: the flexural (ZA) mode and
Kohn anomalies, on the relaxed cell.
