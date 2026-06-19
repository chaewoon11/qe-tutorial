---
title: "Chapter 3 — Brillouin-zone sampling & smearing"
sidebar_label: "3 · BZ sampling & smearing"
sidebar_position: 4
---

# Chapter 3 — Brillouin-zone sampling & smearing

Chapter 2 converged the *basis* (`ecutwfc`/`ecutrho`). The other half of every
DFT calculation's accuracy is how finely you sample the **Brillouin zone** — the
`K_POINTS` line. This chapter explains Monkhorst–Pack grids, why **metals** are
so much harder to converge than semiconductors, and the **smearing** trick that
makes them tractable.

**Prerequisites:** [Chapter 2](./02-pseudopotentials.md).

---

## 1. Why we sum over k-points

Ground-state quantities are integrals over the Brillouin zone. The number of
electrons, the density, the energy — all have the form

$$
\bar{g} = \sum_n \int_\text{BZ} \frac{d\mathbf{k}}{V_\text{BZ}}\, f_{n\mathbf{k}}\,
g_{n\mathbf{k}},
$$

where $f_{n\mathbf{k}}$ is the occupation of band $n$ at $\mathbf{k}$. We can't do
the integral analytically, so we replace it by a weighted **sum over a finite
grid** of k-points. Convergence with respect to that grid is just as essential
as convergence with respect to `ecutwfc`.

## 2. Monkhorst–Pack grids

The standard recipe is a uniform **Monkhorst–Pack** grid: `K_POINTS automatic`
with $N_1 \times N_2 \times N_3$ subdivisions and an optional shift.

```text
K_POINTS automatic
  6 6 6 0 0 0      ! 6×6×6 grid, no shift (Γ-centred)
```

A few practical rules:

- **Match the grid to the cell.** Reciprocal-space spacing scales as $1/a_i$, so
  a long real-space axis needs *fewer* k-points along it. For a cubic cell, use
  an equal $N\times N\times N$.
- **Respect symmetry points.** Choose $N$ so that important high-symmetry points
  fall *on* the grid (e.g. multiples of the relevant integer for hexagonal
  cells). For cubic GaAs an even $N$ includes the zone-boundary X point.
- **Symmetry does the rest.** As you saw in Chapter 1, QE folds the full grid to
  the irreducible wedge, so a $6\times6\times6=216$ grid costs only 16 k-points.

The shift (`0 0 0` vs `1 1 1`) moves the grid off Γ; a shifted grid samples more
inequivalent points for the same $N$, but an unshifted Γ-centred grid is safer
when you later need Γ itself (e.g. for phonons). We use unshifted grids
throughout.

## 3. Metals: the Fermi-surface problem

For an **insulator or semiconductor** (like GaAs), every band is either
completely full ($f=1$) or completely empty ($f=0$) at every k-point. The
integrand is smooth across the whole BZ, so the k-sum converges quickly.

A **metal** is different. A partially-filled band crosses the Fermi level
$\varepsilon_F$, so the occupation

$$
f_{n\mathbf{k}} = \theta(\varepsilon_F - \varepsilon_{n\mathbf{k}})
$$

**jumps from 1 to 0** at the Fermi surface — a discontinuity sitting somewhere
*inside* the BZ. Sampling a step function on a finite grid converges agonizingly
slowly (and erratically): tiny shifts of the grid move points across the Fermi
surface and the energy lurches. With strictly fixed occupations, a metal k-sum
may not converge usefully at all.

## 4. Smearing: soften the step

The cure is to replace the sharp step with a **smooth occupation function** of
width `degauss` (an effective electronic temperature), set with
`occupations = 'smearing'`:

```fortran
&system
    occupations = 'smearing'
    smearing    = 'mv'      ! Marzari–Vanderbilt 'cold' smearing
    degauss     = 0.02      ! width in Ry (~0.27 eV)
/
```

Now $f_{n\mathbf{k}}$ goes smoothly from 1 to 0 over an energy window `degauss`,
so the integrand is smooth and the k-sum converges much faster. The flavours:

| `smearing` | Function | Notes |
|---|---|---|
| `gaussian` (`gauss`) | Gaussian | simplest, large $T S$ correction |
| `mp` | Methfessel–Paxton | small energy error, occupations can go slightly negative |
| `mv` | Marzari–Vanderbilt "cold" | small error, well-behaved — a good default for metals |
| `fd` | Fermi–Dirac | a *physical* temperature $T$ (use when you actually mean finite-$T$) |

`degauss` is a convergence parameter, not a free knob: too small and you lose
the benefit (slow k-convergence again); too large and you smear away real
physics. The correct procedure is to converge the energy jointly in `degauss`
**and** the k-grid (a finer grid tolerates a smaller `degauss`). QE reports a
"smearing contrib. (-TS)" term; the cold/MP schemes keep it tiny so the smeared
energy is close to the true $T\!=\!0$ value.

:::warning Insulator? Use `fixed`
Smearing is for metals. For a semiconductor or insulator with a real gap, use
`occupations = 'fixed'` (as we do for GaAs). Smearing a gapped system just adds
noise and can spuriously partially-occupy the conduction band.
:::

## 5. Hands-on: converge the k-grid for both

Sweep the grid density $N$ for GaAs (insulator, `fixed`) and for fcc aluminium
(metal, `mv` smearing), recording the total energy each time. For GaAs, reuse
the Chapter 0 input and edit the `K_POINTS` line:

```bash
cd code/00-first-scf
for N in 2 4 6 8 10 12; do
  sed -e "s/conv_thr      = 1.0d-20/conv_thr      = 1.0d-10/" \
      -e "s/  6 6 6 0 0 0/  $N $N $N 0 0 0/" inputs/gaas.scf.in > /tmp/k.in
  echo "$N $(grep '^!' <(pw.x < /tmp/k.in) | awk '{print $5}')"
done
```

The aluminium input (a one-atom fcc metal with smearing) is in
[`code/03-brillouin-zone/`](https://github.com/chaewoon11/qe-tutorial/tree/master/code/03-brillouin-zone),
along with both data files.
[`notebooks/03-kpoint-convergence.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/03-kpoint-convergence.ipynb)
plots the two convergence curves together:

![k-point convergence: metal vs. semiconductor](/img/03-kpoint-convergence.png)

The contrast is the lesson of the chapter:

- **GaAs (insulator)** converges fast and smoothly — an $8\times8\times8$ grid
  is already within ~2 meV/atom and $10\times10\times10$ reaches ~0.5 meV/atom,
  because full bands make the BZ integrand smooth.
- **Aluminium (metal)** converges **slowly and non-monotonically** — it is still
  bouncing around the ~5 meV/atom level out to $16\times16\times16$ and only
  settles near $20\times20\times20$, even *with* smearing, because the Fermi
  surface must be finely resolved. A metal routinely wants several times as many
  k-points per direction as an insulator.

:::tip Rule of thumb
Converge the k-grid for the property you report, on the *relaxed* cell, and
re-check it whenever you change the cell size. Metals, small cells, and
fine-energy-scale quantities (magnetism, phonons, transport) demand the densest
grids.
:::

## Exercises

1. **Pick the grids.** From the figure, what $N\times N\times N$ converges each
   material to ~1 meV/atom? How many times denser does the metal need to be?
2. **Smearing matters.** Re-run aluminium at a coarse grid (say $4\times4\times4$)
   with `occupations = 'fixed'` instead of smearing. What happens — does it even
   converge, and how does the energy compare?
3. **`degauss` sensitivity.** At a fixed Al grid, run `degauss` = 0.01, 0.02,
   0.05 Ry. How much does the total energy move? Which way does too-large
   `degauss` bias it?
4. **Cell size and k-density.** If you doubled the GaAs cell along one axis
   (making it twice as long in real space), what should happen to the number of
   k-points you need along that axis, and why?

<details>
<summary>Show solutions</summary>

**1. Pick the grids.** Reading the figure, GaAs reaches ~1 meV/atom by about
$10\times10\times10$ ($8\times8\times8$ already gives ~2 meV/atom), while
aluminium is still oscillating at the ~5 meV/atom level out to
$16\times16\times16$ and only settles near $20\times20\times20$ — several times
denser per direction. The metal's Fermi surface is the reason.

**2. Smearing off on a metal.** With `fixed` occupations and a coarse grid, the
aluminium SCF either fails to converge or lands on a meaningless energy: there is
no gap, so deciding which states are occupied is ambiguous and grid-dependent.
This is exactly the Fermi-surface discontinuity from §3 — smearing exists to fix
it.

**3. `degauss` sensitivity.** The total energy shifts by up to tens of meV across
that range. With Gaussian smearing the energy is pushed *down* (the entropy term
over-stabilises); cold/MP smearing keeps the bias small. The right value is the
largest `degauss` whose energy still agrees with the dense-grid, small-`degauss`
limit — converge the two together.

**4. Cell size and k-density.** Reciprocal-lattice spacing along an axis scales
as $1/a$. Doubling the real-space length halves the reciprocal spacing, so you
need about **half** as many k-points along that axis to keep the same sampling
density. Big cells → coarse grids; small cells → dense grids.

</details>

---

**Next:** Chapter 4 — Structural optimization *(coming soon)*: forces, stress,
and relaxing the cell with `relax` / `vc-relax`.
