---
title: "Chapter 0 — What is Quantum ESPRESSO? Your first SCF"
sidebar_label: "0 · Your first SCF"
sidebar_position: 1
---

# Chapter 0 — What is Quantum ESPRESSO? Your first SCF

**Goal:** understand what a DFT code actually computes, learn the QE file
ecosystem, then run a complete self-consistent-field (SCF) calculation on
silicon and read every important number out of the output.

**Prerequisites:** [Setup](../setup/installation) done (`pw.x` on your `PATH`).

**Time:** ~20 minutes. The calculation itself runs in ~1 second.

---

## 1. The 60-second theory

Density functional theory (DFT) replaces the impossible many-electron
Schrödinger equation with a set of *single-particle* equations — the
**Kohn–Sham (KS) equations** — that have the same ground-state electron density
$n(\mathbf{r})$:

$$
\left[-\tfrac{1}{2}\nabla^2 + v_{\text{KS}}[n](\mathbf{r})\right]\psi_i(\mathbf{r})
= \varepsilon_i\,\psi_i(\mathbf{r}),
\qquad
n(\mathbf{r}) = \sum_i^{\text{occ}} |\psi_i(\mathbf{r})|^2
$$

The catch: the potential $v_{\text{KS}}[n]$ depends on the density $n$, which
depends on the orbitals $\psi_i$, which depend on the potential. So you **guess
a density, build the potential, solve for orbitals, rebuild the density, and
repeat until it stops changing** — that loop is the *self-consistent field*
(SCF) cycle. That is literally what `pw.x` with `calculation='scf'` does.

In a crystal, Bloch's theorem lets us label states by a wavevector $\mathbf{k}$
in the Brillouin zone and expand each orbital in **plane waves** — hence
"plane-wave DFT". Two knobs control accuracy: how many plane waves (`ecutwfc`)
and how finely we sample $\mathbf{k}$ (`K_POINTS`). Those are Chapters 2 and 3;
here we just pick reasonable values and run.

## 2. The QE file ecosystem

| Thing | What it is |
|---|---|
| `*.in` | your **input** — namelists (`&control`, `&system`, …) + cards (`ATOMIC_SPECIES`, …) |
| `*.out` | the **output** — human-readable log; you `grep` results out of it |
| `*.UPF` | a **pseudopotential**, one per element (Unified Pseudopotential Format) |
| `outdir/` (here `./out`) | scratch: wavefunctions, charge density, the `*.save` directory |
| `prefix` | a name tag (`'si'`) prepended to all scratch files so runs don't collide |

You write the `.in`, run `pw.x`, read the `.out`. Everything else is bookkeeping.

## 3. The input file, line by line

Here is [`si.scf.in`](https://github.com/chaewoon11/qe-tutorial/blob/master/code/00-first-scf/inputs/si.scf.in) —
silicon, diamond structure, 2 atoms per cell:

```fortran title="code/00-first-scf/inputs/si.scf.in"
&control
    calculation   = 'scf'          ! self-consistent ground state
    prefix        = 'si'           ! tag for scratch files
    outdir        = './out'        ! where scratch goes
    pseudo_dir    = '../pseudos'   ! where the .UPF lives (relative to THIS dir)
    verbosity     = 'high'         ! print the full energy breakdown
/
&system
    ibrav         = 0              ! lattice given explicitly in CELL_PARAMETERS
    celldm(1)     = 10.26          ! alat: length unit (bohr) for the cell & alat positions
    nat           = 2              ! atoms in the cell
    ntyp          = 1              ! distinct species (just Si)
    ecutwfc       = 40.0           ! plane-wave cutoff (Ry) — the accuracy knob (Ch. 2)
    occupations   = 'fixed'        ! insulator: filled bands, no smearing (Ch. 3)
/
&electrons
    conv_thr      = 1.0d-8         ! stop the SCF loop when energy is this converged (Ry)
    mixing_beta   = 0.7            ! how aggressively to mix new density into old
/
ATOMIC_SPECIES
  Si  28.0855  Si_ONCV_PBE-1.1.upf   ! symbol, mass, pseudopotential file

CELL_PARAMETERS alat                 ! fcc primitive vectors, in units of celldm(1)
  -0.50   0.00   0.50
   0.00   0.50   0.50
  -0.50   0.50   0.00

ATOMIC_POSITIONS alat                ! Cartesian, in units of celldm(1)
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25

K_POINTS automatic                   ! Monkhorst–Pack grid (Ch. 3)
  6 6 6 0 0 0
```

:::tip Why `ibrav = 0`?
You *could* let QE build the fcc lattice from a symmetry code (`ibrav = 2`), but
writing `CELL_PARAMETERS` explicitly is more transparent and generalizes to any
structure — so we use it throughout this tutorial.
:::

## 4. Run it

From the `code/00-first-scf/` directory (so the relative `pseudo_dir` resolves):

```bash
cd code/00-first-scf
pw.x < inputs/si.scf.in > my.scf.out
```

That's it. A reference copy of the output is in
[`ref/si.scf.out`](https://github.com/chaewoon11/qe-tutorial/blob/master/code/00-first-scf/ref/si.scf.out) —
compare yours against it.

:::note Serial vs. parallel
`pw.x < in > out` runs serial. To use 4 MPI ranks:
`mpirun -np 4 pw.x < inputs/si.scf.in > my.scf.out`. Parallelization is
Chapter 12 — for a 2-atom cell, serial is plenty.
:::

## 5. Read the output

Don't be intimidated by 1000 lines — you `grep` for the few that matter.

**Did it finish and converge?**

```bash
grep -E "JOB DONE|convergence has been achieved" my.scf.out
```

```text
     convergence has been achieved in   6 iterations
   JOB DONE.
```

The SCF loop reached `conv_thr` in 6 cycles. (No "JOB DONE" → it crashed; read
the bottom of the file for the error.)

**Watch the SCF loop converge** — this is the self-consistency from §1 in action:

```bash
grep "estimated scf accuracy" my.scf.out
```

```text
     estimated scf accuracy    <       0.05325399 Ry
     estimated scf accuracy    <       0.00243937 Ry
     estimated scf accuracy    <       0.00004851 Ry
     estimated scf accuracy    <       0.00000053 Ry
     estimated scf accuracy    <       0.00000002 Ry     ← below conv_thr (1e-8) → done
```

Each iteration the density gets closer to self-consistent and the error drops
by ~1–2 orders of magnitude. That funnel *is* the algorithm working.

![Silicon SCF convergence](/img/00-scf-convergence.png)

:::tip Companion notebook
The plot above is produced by
[`notebooks/00-scf-convergence.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/00-scf-convergence.ipynb),
which parses the reference output and plots the accuracy vs. iteration. Grab it
to make the same figure from *your* run.
:::

**The headline number — total energy:**

```bash
grep '^!' my.scf.out
```

```text
!    total energy              =     -15.76393891 Ry
```

The `!` marks the final, converged total energy. With `verbosity='high'` you
also get its breakdown:

```text
     one-electron contribution =       4.75014667 Ry
     hartree contribution      =       1.10501663 Ry   (classical e–e repulsion)
     xc contribution           =      -4.81817260 Ry   (exchange–correlation)
     ewald contribution        =     -16.80092960 Ry   (ion–ion repulsion)
```

:::warning Absolute energies are meaningless on their own
A total energy depends on the pseudopotential. Only *differences* (between
structures, volumes, phases) are physical. We exploit that constantly from
Chapter 3 on.
:::

**A physical observable — the highest occupied level:**

```bash
grep "highest occupied" my.scf.out
```

```text
     highest occupied level (ev):     6.2395
```

With 8 valence electrons filling 4 bands, this is the valence-band maximum. (In
Chapter 5 we turn this into a real band structure and read off the gap.)

## 6. What just happened

You guessed a density, QE ran the Kohn–Sham SCF loop until the density stopped
changing, and reported the ground-state energy of crystalline silicon. Every
later calculation — relaxations, bands, phonons — is built on this exact loop.

## Exercises

1. **Convergence behaviour.** Re-run with `conv_thr = 1.0d-4`. How many SCF
   iterations now? What's the total energy, and how does it differ from the
   `1.0d-8` result?
2. **The mixing knob.** Set `mixing_beta = 0.3`. Does the SCF take more or fewer
   iterations? Why might a *smaller* mixing be safer for harder systems?
3. **Break it on purpose.** Change the pseudopotential filename to something
   that doesn't exist and run. Read the error message — get used to what a
   failure looks like.
4. **Read the geometry.** Find the lines reporting `unit-cell volume` and
   `number of Kohn-Sham states`. Why are there exactly 4 KS states here?

<details>
<summary>Show solutions</summary>

**1. Convergence behaviour (`conv_thr = 1.0d-4`)**

- **3 iterations** (vs. 6 with `1.0d-8`) — a looser target is reached sooner.
- Total energy: **−15.76389664 Ry**, vs. **−15.76393891 Ry** at `1.0d-8`.
- Difference: **~4×10⁻⁵ Ry ≈ 0.6 meV**. Tiny here, but the total energy is only
  as converged as `conv_thr` lets it be. For *energy differences* you want both
  runs converged well below the difference you care about — `1.0d-6` Ry or
  tighter is the usual safe default.

**2. The mixing knob (`mixing_beta = 0.3`)**

- **14 iterations** (vs. 6 at `beta = 0.7`) — same final energy
  (**−15.76393888 Ry**), just slower.
- `mixing_beta` controls how much of the *new* density is blended into the old
  one each cycle. Large beta = bigger steps = faster *when it works*, but can
  oscillate or diverge for "hard" systems (metals, magnetic materials, big
  cells). A **smaller beta takes more, smaller, safer steps** — the standard fix
  when an SCF won't converge. Silicon is easy, so 0.7 is fine and faster.

**3. Break it on purpose**

With a bad pseudopotential filename, `pw.x` stops almost immediately:

```text
     Error in routine readpp (1):
     file .../Si_does_not_exist.upf not found

     stopping ...
```

QE prints errors in a `%%%%`-bordered block naming the **routine** (`readpp`)
and the reason, then `stopping ...`. There is **no "JOB DONE"**. First debugging
move on any failed run: read the last ~20 lines of the output.

**4. Read the geometry**

```bash
grep -E "unit-cell volume|number of Kohn-Sham states" my.scf.out
```

```text
     unit-cell volume          =     270.0114 (a.u.)^3
     number of Kohn-Sham states=            4
```

Silicon contributes **4 valence electrons per atom × 2 atoms = 8 electrons**.
This is a (non-magnetic) insulator with doubly-occupied states, so it needs
8 / 2 = **4 occupied bands** — hence 4 KS states by default. (When we add empty
bands for a band structure in Chapter 5, we raise this with `nbnd`.)

</details>

---

**Next:** Chapter 1 — Plane-wave DFT & the SCF loop *(coming soon)*, where we
open up the theory behind every line you just ran.
