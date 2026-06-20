---
id: 2d-cell
title: "A1 — Building a 2D cell (graphene)"
sidebar_label: "A1 · Building a 2D cell"
sidebar_position: 1
---

# A1 — Building a 2D cell (graphene)

The core tutorial used 3D crystals (GaAs, Fe, NiO). This advanced track turns to
**2D materials**, with **graphene** as the running example. A plane-wave code is
inherently 3D-periodic, so the first job is to *fake* a single isolated sheet
correctly: a slab in a box with vacuum, plus a boundary correction. We also meet
the **2D Brillouin zone** and why graphene needs smearing.

**Prerequisites:** the core tutorial — especially [BZ sampling](../chapters/brillouin-zone.md)
and [SCF convergence](../chapters/first-scf.md).

---

:::note The workflow — do it the proper way
Every property from here on is computed **on the same converged, relaxed
structure**, not on a hand-typed guess. The discipline:

1. **Converge the parameters.** For a 2D cell the new one is the **vacuum**; we
   fix a reasonable `ecutwfc/ecutrho` and a `12×12×1` k-grid and test only the
   vacuum (§3).
2. **Relax** the structure with those fixed parameters (§4).
3. **Reuse** that relaxed cell for *all* downstream calculations (bands, phonons…).

And we use a consistent set of self-consistency thresholds throughout the 2D
track:

| Calculation | threshold | value |
|---|---|---|
| `scf` | `conv_thr` | `1.0d-20` |
| `relax` / `vc-relax` | `conv_thr` | `1.0d-18` |
| `nscf` (bands/DOS) | `conv_thr` | `1.0d-16` |
| `ph.x` (phonons) | `tr2_ph` | `1.0d-20` |
:::

## 1. A 2D sheet in a 3D code: vacuum + a boundary fix

QE imposes periodic boundary conditions in all three directions. To model **one**
graphene sheet we use a **supercell**: the atoms lie in the $xy$-plane and we add
**vacuum** along $z$ so the sheet doesn't see its periodic images above and below.

Vacuum alone has a subtle problem: even far apart, the periodic images interact
**electrostatically** (a charged or polarized slab sets up a long-range field that
the periodic boundary conditions wrap around). For a neutral, non-polar sheet like
graphene this is weak, but the clean fix is **`assume_isolated = '2D'`**, which
applies a *cutoff Coulomb* in the out-of-plane direction so the sheet is truly
isolated. Two requirements:

```fortran
&system
    assume_isolated = '2D'   ! cutoff Coulomb along z (slab must be xy, centred at z=0)
    ...
/
```

- the slab normal is the **third cell axis** ($z$), and
- the atoms are **centred at $z = 0$** (so the vacuum is split symmetrically).

We build graphene with an explicit hexagonal cell (keeping the `ibrav = 0`
convention of the core tutorial), the two carbon atoms at $z=0$:

Here is the complete canonical input — the **relaxed** cell (§4) with the
converged vacuum (§3) and the `scf` threshold from the convention above:

```fortran title="code/advanced/01-2d-cell/graphene.scf.in"
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
    assume_isolated = '2D'        ! cutoff Coulomb along z (slab in xy, centred at z=0)
    occupations = 'smearing'      ! graphene is a semimetal
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
   2.467137   0.000000   0.000000     ! a (relaxed, §4)
  -1.233568   2.136488   0.000000     ! a2 = a(-1/2, √3/2, 0)
   0.000000   0.000000  20.000000     ! c = vacuum thickness (§3)
ATOMIC_POSITIONS crystal
  C 0.000000000 0.000000000 0.0       ! both atoms at z = 0
  C 0.333333333 0.666666667 0.0
K_POINTS automatic
  12 12 1 0 0 0                        ! N×N×1, N a multiple of 3 (K on grid)
```

(The vacuum test in §3 starts from the experimental $a=2.46$ Å; §4 then relaxes
$a$ to 2.4671 Å. Vacuum convergence is insensitive to that ~0.3 % change.)

## 2. The 2D Brillouin zone and smearing

Graphene's BZ is a **hexagon**; the standard path is **Γ–M–K–Γ**, and the physics
lives at the corner **K**. Two practical points:

- **k-grid `N × N × 1`.** Only sample in-plane — one k-point along the
  out-of-plane direction (there is no dispersion through vacuum). For graphene the
  important **K** point $(\tfrac13,\tfrac13,0)$ falls *on* a Γ-centred grid only
  when **`N` is a multiple of 3**; we use `12 12 1`.
- **Smearing.** Graphene is a **semimetal** — valence and conduction bands touch
  at K (the Dirac point), so there is no gap and you must use
  `occupations = 'smearing'` (as for a metal, Chapter 3), e.g. `mv` with a small
  `degauss`.

## 3. Hands-on: converge the vacuum

The one new convergence parameter for a 2D cell is the **vacuum thickness** $c$.
Too little and the sheet feels its images; too much wastes plane waves (a taller
box at fixed `ecutrho` means more G-vectors). Sweep it and watch the total energy
settle:

```bash
cd code/advanced/01-2d-cell
for c in 10 14 18 22 26 30; do
  sed "s/__C__/$c.0/" graphene.scf.tmpl > /tmp/g.in   # set the c lattice vector
  echo "$c $(grep '^!' <(pw.x < /tmp/g.in) | tail -1 | awk '{print $5}')"
done
```

[`notebooks/A1-vacuum.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/A1-vacuum.ipynb)
plots the result:

![Graphene vacuum convergence](/img/A1-vacuum.png)

With `assume_isolated='2D'` the total energy is flat to **under 0.03 meV/atom**
across 14–30 Å — the sheet is effectively isolated already at ~14–15 Å (without
the flag you would need far more vacuum). But there is a **hard floor**: a box
of **`c ≲ 10 Å` errors out** with `2D cutoff smaller than pseudo cutoff radius`,
because the cell must be taller than the pseudopotential's cutoff radius. We adopt
**`c = 20 Å`** — safely above the floor and fully converged.

:::tip `assume_isolated='2D'` lets you use *less* vacuum
With the cutoff-Coulomb correction the energy converges at a smaller $c$ than
plain vacuum would need, because the spurious image interaction is removed
analytically rather than by brute-force separation. Still converge it — the
correction is exact for electrostatics, not for wavefunction overlap.
:::

## 4. Relax, then reuse

With the parameters fixed (`ecutwfc=50`, `ecutrho=400`, `12×12×1`, `c=20 Å`), we
relax the in-plane lattice. Graphene's two atoms are fixed by symmetry, so only
the lattice constant can change — a **`vc-relax` with `cell_dofree='2Dxy'`**
(relax $a$, **keep the vacuum $c$ fixed**), at the relax threshold `conv_thr=1.0d-18`:

```fortran title="code/advanced/01-2d-cell/graphene.relax.in"
&control
    calculation = 'vc-relax'
    prefix      = 'graphene'
    outdir      = './out'
    pseudo_dir  = '../../pseudos'
    verbosity   = 'high'
    tprnfor = .true.
    tstress = .true.
    forc_conv_thr = 1.0d-4
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
    conv_thr = 1.0d-18
    mixing_beta = 0.7
/
&ions
    ion_dynamics = 'bfgs'
/
&cell
    cell_dofree = '2Dxy'          ! relax in-plane lattice only; keep vacuum c fixed
    press = 0.0d0
/
ATOMIC_SPECIES
  C 12.011 C.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
   2.460000   0.000000   0.000000     ! start from the experimental a
  -1.230000   2.130000   0.000000
   0.000000   0.000000  20.000000
ATOMIC_POSITIONS crystal
  C 0.000000000 0.000000000 0.0
  C 0.333333333 0.666666667 0.0
K_POINTS automatic
  12 12 1 0 0 0
```

```bash
cd code/advanced/01-2d-cell
pw.x < graphene.relax.in > graphene.relax.out
grep -A4 "Begin final coordinates" graphene.relax.out
```

It converges in 3 BFGS steps to a relaxed lattice constant

$$
a_0^{\text{PBE}} = 2.4671\ \text{Å} \quad (\text{C–C bond } a_0/\sqrt3 = 1.424\ \text{Å}),
$$

versus the experimental **2.46 Å** — PBE's usual slight over-estimate (~0.3 %).
**This relaxed cell** (written into `graphene.scf.in`) is the structure every later
chapter — Dirac bands (A2), phonons (A4) — starts from. Relax once, reuse
everywhere; never recompute a property on an unrelaxed guess.

## Exercises

1. **Why `N` multiple of 3?** Sketch the hexagonal BZ and a Γ-centred 12×12 grid.
   Show that K $(\tfrac13,\tfrac13)$ is a grid point but that an 8×8 grid misses
   it. Why does missing K matter for graphene specifically?
2. **Vacuum vs cost.** Doubling $c$ at fixed `ecutrho` roughly doubles the number
   of G-vectors. From your sweep, what is the smallest $c$ within ~1 meV/atom of
   the converged energy — and why is "more vacuum" not free?
3. **Drop the 2D flag.** Re-run at small $c$ with and without
   `assume_isolated='2D'`. Which converges to the right answer with less vacuum,
   and what is the flag physically correcting?
4. **Center matters.** Move both atoms to $z = 0.5$ (crystal) instead of $0$ and
   re-run with `assume_isolated='2D'`. What goes wrong, and why does the slab need
   to sit at $z=0$?

<details>
<summary>Show solutions</summary>

**1. Multiple of 3.** A Γ-centred $N\times N$ mesh has points at $(i/N, j/N)$.
$K=(1/3,1/3)$ is on the mesh iff $N/3$ is an integer. 8×8 has no point at exactly
$1/3$, so it misses the Dirac point — bad for graphene because *all* the
low-energy physics (the Dirac cone, the density of states at $E_F$) is at K.

**2. Vacuum vs cost.** Read the smallest $c$ whose energy is within ~1 meV/atom of
the $c=30$ Å value (typically ~15–18 Å with `assume_isolated='2D'`). More vacuum
isn't free because the FFT box grows with $c$: at fixed `ecutrho` the number of
$z$ G-vectors scales with the box height, so the calculation gets slower and uses
more memory for empty space.

**3. Drop the 2D flag.** With `assume_isolated='2D'` the energy converges at
smaller $c$; plain vacuum needs more separation to suppress the periodic-image
electrostatics. The flag replaces that brute-force separation with an exact
cutoff-Coulomb boundary condition along $z$.

**4. Center matters.** `assume_isolated='2D'` assumes the slab is centred at
$z=0$; putting the atoms at $z=0.5$ places them at the cell edge, splitting the
sheet across the periodic boundary and putting the vacuum (and the Coulomb
cutoff) in the wrong place — the correction is applied through the middle of the
sheet. Always centre a 2D slab at $z=0$.

</details>

---

**Next:** [A2 — Dirac bands of graphene](./02-dirac-bands.md): the band structure
along Γ–M–K–Γ and the Dirac cone at K.
