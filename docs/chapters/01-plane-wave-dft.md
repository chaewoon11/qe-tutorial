---
title: "Chapter 1 — Plane-wave DFT & the SCF loop"
sidebar_label: "1 · Plane-wave DFT & SCF"
sidebar_position: 2
---

# Chapter 1 — Plane-wave DFT & the SCF loop

In Chapter 0 you ran an SCF and read off a total energy. Now we open the box:
*what* is being solved, *why* in plane waves, and *how* the self-consistent loop
actually works. Every concept here maps to a line you already saw in the GaAs
output — we'll keep pointing back to it.

**Prerequisites:** [Chapter 0](./00-first-scf.md) (you have a `gaas.scf.out` to look at).

---

## 1. From the many-body problem to Kohn–Sham

A solid is a horde of nuclei and electrons. Fixing the nuclei
(Born–Oppenheimer), the electronic Hamiltonian is

$$
\hat{H} = \underbrace{-\tfrac{1}{2}\sum_i \nabla_i^2}_{\text{kinetic}}
        + \underbrace{\sum_i v_\text{ext}(\mathbf{r}_i)}_{\text{nuclei}\to\text{electrons}}
        + \underbrace{\tfrac{1}{2}\sum_{i\neq j}\frac{1}{|\mathbf{r}_i-\mathbf{r}_j|}}_{\text{electron–electron}} .
$$

That last term couples all electrons, so the $3N$-dimensional wavefunction
$\Psi(\mathbf{r}_1,\dots,\mathbf{r}_N)$ is hopeless to store, let alone solve, for
a real crystal.

**Hohenberg–Kohn** rescue: the ground-state energy is a *functional of the
electron density* $n(\mathbf{r})$ alone — a function of just 3 variables — and it
is minimised by the true density. **Kohn–Sham** make it practical: replace the
interacting system with a fictitious system of *non-interacting* electrons that
has the **same density**. Those electrons obey the single-particle equations
from Chapter 0,

$$
\left[-\tfrac{1}{2}\nabla^2 + v_\text{KS}(\mathbf{r})\right]\psi_i = \varepsilon_i\,\psi_i,
\qquad
v_\text{KS} = v_\text{ext} + \underbrace{v_\text{H}[n]}_{\text{Hartree}} + \underbrace{v_\text{xc}[n]}_{\text{exchange–correlation}} .
$$

The three pieces are explicit functionals of the geometry and density. The
**external potential** is the electrostatic attraction of the nuclei (charge
$Z_I$ at position $\mathbf{R}_I$) acting on the electrons,

$$
v_\text{ext}(\mathbf{r}) = -\sum_I \frac{Z_I}{|\mathbf{r}-\mathbf{R}_I|},
$$

which in practice is replaced by the smooth **pseudopotential** for each species
(§5). The **Hartree** term is the classical electrostatic potential of the
electron density itself,

$$
v_\text{H}(\mathbf{r}) = \int \frac{n(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}\,d\mathbf{r}',
$$

and the **exchange–correlation** potential is the functional derivative
$v_\text{xc}(\mathbf{r}) = \delta E_\text{xc}[n]/\delta n(\mathbf{r})$ — the one
piece we must approximate. Only $v_\text{ext}$ knows what material this is;
$v_\text{H}$ and $v_\text{xc}$ depend on $\mathbf{r}$ only through $n$, which is
exactly why the equations must be solved self-consistently.

Correspondingly, the total energy splits into terms you have already met in the
output:

```text
     one-electron contribution =    -120.23947851 Ry   ! kinetic + v_ext on the electrons
     hartree contribution      =      45.04959364 Ry   ! v_H : classical e–e repulsion
     xc contribution           =     -50.87806355 Ry   ! E_xc : everything quantum-mechanical
     ewald contribution        =     -93.01921000 Ry   ! nucleus–nucleus (Ewald sum)
     one-center paw contrib.   =    -234.87770454 Ry   ! PAW augmentation inside the cores
```

Everything hard about the interacting electron gas is swept into the single
**exchange–correlation** term $E_\text{xc}[n]$, which we don't know exactly and
must approximate (LDA, PBE, … — Chapter 8). The `Functional: SLA PW PBX PBC`
line near the top of `gaas.scf.out` records that this run used **PBE**.

## 2. Bloch's theorem and the plane-wave basis

In a periodic crystal the potential satisfies $v(\mathbf{r}+\mathbf{R}) =
v(\mathbf{r})$ for any lattice vector $\mathbf{R}$. **Bloch's theorem** then says
every eigenstate can be labelled by a wavevector $\mathbf{k}$ and written as a
plane wave times a cell-periodic function:

$$
\psi_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}}\, u_{n\mathbf{k}}(\mathbf{r}),
\qquad u_{n\mathbf{k}}(\mathbf{r}+\mathbf{R}) = u_{n\mathbf{k}}(\mathbf{r}).
$$

Because $u_{n\mathbf{k}}$ is periodic, it expands exactly in **plane waves**
whose wavevectors are reciprocal-lattice vectors $\mathbf{G}$:

$$
\psi_{n\mathbf{k}}(\mathbf{r}) = \sum_{\mathbf{G}} c_{n\mathbf{k}}(\mathbf{G})\,
e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}} .
$$

This is the entire reason QE is a "plane-wave" code. The infinite sum is
truncated by keeping only plane waves with kinetic energy below the cutoff —
this is exactly what `ecutwfc` sets:

$$
\tfrac{1}{2}\,|\mathbf{k}+\mathbf{G}|^2 < E_\text{cut} \;(=\texttt{ecutwfc}).
$$

So the basis is a sphere in $\mathbf{G}$-space of radius $\sqrt{2E_\text{cut}}$.
The number of plane waves is the number of $\mathbf{G}$ inside it, which scales
as the sphere's volume:

$$
N_\text{PW} \;\propto\; E_\text{cut}^{3/2}.
$$

### Hands-on: watch the basis grow

Run the Chapter 0 input at several cutoffs and grep the basis size:

```bash
cd code/00-first-scf
for ec in 20 30 40 50 60 70 80; do
  sed "s/ecutwfc       = 50.0/ecutwfc       = ${ec}.0/" inputs/gaas.scf.in > /tmp/s.in
  pw.x < /tmp/s.in > /tmp/s.out
  echo "$ec $(grep 'Smooth grid' /tmp/s.out | grep -oE '[0-9]+ G-vectors')"
done
```

The collected numbers are in
[`code/01-plane-wave-dft/npw_sweep.dat`](https://github.com/chaewoon11/qe-tutorial/blob/master/code/01-plane-wave-dft/npw_sweep.dat),
and [`notebooks/01-basis-size.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/01-basis-size.ipynb)
plots them:

![Plane-wave basis size vs. cutoff](/img/01-basis-size.png)

The points lie exactly on the $E_\text{cut}^{3/2}$ line: from 20 → 80 Ry,
`ecutwfc` ×4 makes the basis ×8 (3695 → 29603 plane waves).

:::warning Basis *size* ≠ accuracy
A bigger basis is more expensive but not automatically "better converged."
Choosing `ecutwfc` so that the *energy* (or force, or gap) is converged is a
separate study — Chapter 3.
:::

### Two grids: smooth and dense

Look near the top of `gaas.scf.out`:

```text
     Dense  grid:    41151 G-vectors     FFT dimensions: (  50,  50,  50)
     Smooth grid:    14361 G-vectors     FFT dimensions: (  36,  36,  36)
```

The **smooth grid** (set by `ecutwfc`) holds the *wavefunctions*; the **dense
grid** (set by `ecutrho`) holds the *charge density* $n = \sum|\psi|^2$, which
contains shorter-wavelength features and so needs a finer grid. For
norm-conserving pseudopotentials `ecutrho = 4×ecutwfc` suffices and the two
nearly coincide; for the **PAW** potentials we use here, the density is sharper
and we set `ecutrho` higher (we used 400 = 8×50). QE moves between real and
reciprocal space with **FFTs** — the `FFT dimensions` are those grids.

## 3. Sampling the Brillouin zone

Ground-state quantities are integrals over the Brillouin zone, e.g. the density

$$
n(\mathbf{r}) = \sum_n \int_\text{BZ} \frac{d\mathbf{k}}{V_\text{BZ}}\,
f_{n\mathbf{k}}\,|\psi_{n\mathbf{k}}(\mathbf{r})|^2 ,
$$

approximated by a **finite sum** over a grid of $\mathbf{k}$-points — that's the
`K_POINTS automatic / 6 6 6` line. Crystal symmetry means many $\mathbf{k}$ are
equivalent, so QE only computes the **irreducible** set:

```text
     24 Sym. Ops. (no inversion) found
     ...
     number of k points=    16
```

Two things to notice. First, the $6\times6\times6 = 216$ grid collapses to just
**16 irreducible** points. Second, GaAs has **24** symmetry operations and *no
inversion* — silicon's diamond lattice has 48 *with* inversion. The difference
is physical: in zincblende the two sublattice sites hold *different* atoms (Ga
vs As), which breaks inversion symmetry. (That broken symmetry is also what
makes GaAs piezoelectric and gives it LO–TO phonon splitting — Chapter 11.)
How dense a grid you actually need is Chapter 3.

## 4. The SCF loop, step by step

Now the loop itself. The KS potential depends on the density it is supposed to
produce, so we iterate:

```text
  guess n(r)                         ← superposition of atomic densities
      │
      ▼
  build v_KS[n] = v_ext + v_H[n] + v_xc[n]
      │
      ▼
  solve [ -½∇² + v_KS ] ψ_i = ε_i ψ_i   ← iterative diagonalization (Davidson)
      │                                   for every irreducible k-point
      ▼
  new density  n_out(r) = Σ f_i |ψ_i|²
      │
      ▼
  mix:  n_in ← (1-β) n_in + β n_out      ← mixing_beta = β
      │
      ▼
  converged?  |E_n - E_{n-1}| < conv_thr ?
      │  no → loop back with the mixed n_in
      │  yes ↓
   DONE → total energy, eigenvalues
```

Each box is a Chapter-0 input flag:

| Step | Controlled by | Chapter-0 value |
|---|---|---|
| iterative diagonalization | `diagonalization` | `david` (default) |
| density mixing | `mixing_beta` | `0.7` |
| stop condition | `conv_thr` | `1.0d-20` Ry |

The `estimated scf accuracy` you plotted in Chapter 0 is the loop's own estimate
of how far the energy still is from self-consistency; when it drops below
`conv_thr`, QE stops. Exercise 2 of Chapter 0 (smaller `mixing_beta` → more
iterations) is just this loop taking smaller, safer steps.

:::tip Production rule of thumb
QE's default `conv_thr` (`1.0d-6`) is fine for a quick total energy, but in my
experience a target of about **`1.0d-20` Ry** is a good habit for
production-quality work. It costs only a handful of extra iterations on an easy
system like this, and it pays off later: forces, stresses, and especially phonon
(DFPT) calculations differentiate the energy and are far more sensitive to an
under-converged density than the energy itself is. Converge tight from the start
so the same setup carries through the whole workflow.
:::

### What the solver produces

After "End of self-consistent calculation", QE prints the Kohn–Sham eigenvalues
$\varepsilon_{n\mathbf{k}}$. At $\Gamma$ (the first k-point):

```text
          k = 0.0000 0.0000 0.0000 (  1837 PWs)   bands (ev):

    -7.1922  -7.1922  -7.1922  -7.1139  -7.1139  -5.1109   7.6256   7.6256
     7.6256
```

Read it: **`1837 PWs`** is the actual plane-wave basis for this $\mathbf{k}$
(the per-k count, smaller than the global grid). The **9 eigenvalues** are the 9
bands (= 18 electrons / 2). The lowest five (~−7.1 eV) are the flat **Ga 3d**
semicore bands; then an s-like state; the top **three are degenerate at 7.6256
eV** — the valence-band maximum, whose triple degeneracy is forced by the cubic
symmetry at $\Gamma$. Turning these eigenvalues along a k-path into a band
structure is Chapter 5.

## 5. Why this needs pseudopotentials

A plane-wave basis struggles with the sharp, oscillating wavefunctions near a
nucleus — resolving Ga's $1s$ core would need a cutoff thousands of Ry. The fix
is to replace the core electrons and the strong nuclear potential with a smooth
**pseudopotential**, so the valence wavefunctions are smooth and a ~50 Ry cutoff
works. That is the `*.upf` files in your input, and the `one-center paw contrib.`
energy term is the PAW scheme's bookkeeping for what was smoothed out. Which
pseudopotential to trust, and how it fixes `ecutwfc`/`ecutrho`, is Chapter 2.

## Exercises

1. **Count by hand.** At $\Gamma$, $\mathbf{k}=0$, so the basis is all
   $\mathbf{G}$ with $\tfrac12|\mathbf{G}|^2 < 50$ Ry. Using the cell volume
   from the output, estimate $N_\text{PW}$ and compare with the `1837 PWs` QE
   reports. (Hint: with $E_\text{cut}$ in Ry and $V_\text{cell}$ in bohr³,
   $N_\text{PW} \approx \frac{1}{6\pi^2}\,E_\text{cut}^{3/2}\,V_\text{cell}$.)
2. **Symmetry off.** Add `nosym = .true.` to `&system` and re-run. How many
   k-points now? Does the total energy change? Explain.
3. **Read the bands.** From the $\Gamma$ eigenvalues, which states are the Ga 3d
   semicore manifold and which are the sp valence bands? How many of each, and
   why?
4. **Density needs a finer grid.** Lower `ecutrho` to `200` and re-run. Does the
   `Dense grid` size change? Does the energy? Why is the density grid allowed to
   be coarser or finer independently of `ecutwfc`?

<details>
<summary>Show solutions</summary>

**1. Count by hand.** With $V_\text{cell} = 304.79\ \text{bohr}^3$ and
$E_\text{cut} = 50$ Ry, the plane waves are the reciprocal-lattice points inside
a sphere of radius $\sqrt{E_\text{cut}}$ (in bohr⁻¹, since in Rydberg units the
kinetic energy of $e^{i\mathbf{G}\cdot\mathbf{r}}$ is $|\mathbf{G}|^2$). Their
count is the sphere volume divided by the reciprocal cell volume
$(2\pi)^3/V_\text{cell}$:
$N_\text{PW} \approx \frac{1}{6\pi^2}\,E_\text{cut}^{3/2}\,V_\text{cell}
= \frac{1}{59.2}\times 353.6 \times 304.79 \approx 1.8\times10^3$ — right on
QE's **1837**. The estimate works because the basis is just lattice points in a
sphere.

**2. Symmetry off (`nosym = .true.`).** QE can no longer fold the grid, so it
uses all **216** k-points instead of 16. The total energy is **unchanged** (to
numerical precision): symmetry-equivalent k-points carry identical information,
so summing the full grid gives the same density and energy — just ~13× more
slowly. Symmetry is a pure speed-up here.

**3. Read the bands.** The five eigenvalues near −7.1 eV are the **Ga 3d**
semicore bands (10 electrons → 5 doubly-occupied bands), flat because the 3d
orbitals are tightly bound and barely disperse. The remaining four bands are the
**sp valence** manifold (8 electrons → 4 bands): one s-like state (~−5.1 eV) and
the three degenerate p-like states at the VBM. 5 + 4 = 9 bands = 18/2.

**4. Density grid.** Lowering `ecutrho` to 200 shrinks the `Dense grid` and
changes the energy slightly (the density is now under-resolved — for these PAW
potentials 200 is too low). `ecutrho` is independent of `ecutwfc` because the
wavefunctions and the density they build are *different objects*: the density
$n=\sum|\psi|^2$ has Fourier components up to twice the wavefunction's, and for
PAW/ultrasoft potentials even higher, so it legitimately needs its own,
finer grid.

</details>

---

**Next:** [Chapter 2 — Pseudopotentials & cutoffs](./02-pseudopotentials.md):
what NC / US / PAW actually mean, and how to choose `ecutwfc` and `ecutrho` for
a given set.
