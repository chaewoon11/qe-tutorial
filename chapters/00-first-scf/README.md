# Chapter 0 — What is Quantum ESPRESSO? Your first SCF

**Goal:** understand what a DFT code actually computes, learn the QE file
ecosystem, then run a complete self-consistent-field (SCF) calculation on
silicon and read every important number out of the output.

**Prerequisites:** [`../../setup/`](../../setup/) done (`pw.x` on your `PATH`).

**Time:** ~20 minutes. The calculation itself runs in ~1 second.

---

## 1. The 60-second theory

Density functional theory (DFT) replaces the impossible many-electron
Schrödinger equation with a set of *single-particle* equations — the
**Kohn–Sham (KS) equations** — that have the same ground-state electron density
n(**r**):

$$\left[-\tfrac{1}{2}\nabla^2 + v_{\text{KS}}[n](\mathbf{r})\right]\psi_i = \varepsilon_i\,\psi_i, \qquad n(\mathbf{r}) = \sum_i^{\text{occ}} |\psi_i(\mathbf{r})|^2$$

The catch: the potential $v_{\text{KS}}[n]$ depends on the density $n$, which
depends on the orbitals $\psi_i$, which depend on the potential. So you **guess
a density, build the potential, solve for orbitals, rebuild the density, and
repeat until it stops changing** — that loop is the *self-consistent field*
(SCF) cycle. That is literally what `pw.x calculation='scf'` does.

In a crystal, Bloch's theorem lets us label states by a wavevector **k** in the
Brillouin zone and expand each orbital in **plane waves** — hence "plane-wave
DFT". Two knobs control accuracy: how many plane waves (`ecutwfc`) and how
finely we sample **k** (`K_POINTS`). Those are Chapters 2 and 3; here we just
pick reasonable values and run.

Full derivation: [`../../theory-appendix/`](../../theory-appendix/) *(coming with
Chapter 1)*.

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

Open [`inputs/si.scf.in`](inputs/si.scf.in). Silicon, diamond structure, 2 atoms
per cell:

```fortran
&control
    calculation   = 'scf'          ! self-consistent ground state
    prefix        = 'si'           ! tag for scratch files
    outdir        = './out'        ! where scratch goes
    pseudo_dir    = '../../setup/pseudos'   ! where the .UPF lives (relative to THIS chapter dir)
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

Why `ibrav = 0`? You *could* let QE build the fcc lattice from a symmetry code
(`ibrav = 2`), but writing `CELL_PARAMETERS` explicitly is more transparent and
generalizes to any structure — so we use it throughout this tutorial.

## 4. Run it

From **this chapter's directory** (so the relative `pseudo_dir` resolves):

```bash
cd chapters/00-first-scf
pw.x < inputs/si.scf.in > my.scf.out
```

That's it. A reference copy of the output is in
[`ref/si.scf.out`](ref/si.scf.out) — compare yours against it.

> **Serial vs. parallel.** `pw.x < in > out` runs serial. To use 4 MPI ranks:
> `mpirun -np 4 pw.x < inputs/si.scf.in > my.scf.out`. Parallelization is
> Chapter 12 — for a 2-atom cell, serial is plenty.

## 5. Read the output

Don't be intimidated by 1000 lines — you `grep` for the few that matter.

**Did it finish and converge?**

```bash
grep -E "JOB DONE|convergence has been achieved" my.scf.out
```
```
     convergence has been achieved in   6 iterations
   JOB DONE.
```
The SCF loop reached `conv_thr` in 6 cycles. (No "JOB DONE" → it crashed; read
the bottom of the file for the error.)

**Watch the SCF loop converge** — this is the self-consistency from §1 in action:

```bash
grep "estimated scf accuracy" my.scf.out
```
```
     estimated scf accuracy    <       0.05325399 Ry
     estimated scf accuracy    <       0.00243937 Ry
     estimated scf accuracy    <       0.00004851 Ry
     estimated scf accuracy    <       0.00000053 Ry
     estimated scf accuracy    <       0.00000002 Ry     ← below conv_thr (1e-8 Ry) → done
```
Each iteration the density is closer to self-consistent and the error drops by
~1–2 orders of magnitude. That funnel *is* the algorithm working.

**The headline number — total energy:**

```bash
grep '^!' my.scf.out
```
```
!    total energy              =     -15.76393891 Ry
```
The `!` marks the final, converged total energy. With `verbosity='high'` you
also get its breakdown:
```
     one-electron contribution =       4.75014667 Ry
     hartree contribution      =       1.10501663 Ry   (classical e–e repulsion)
     xc contribution           =      -4.81817260 Ry   (exchange-correlation)
     ewald contribution        =     -16.80092960 Ry   (ion–ion repulsion)
```

> **Absolute total energies are meaningless on their own** — they depend on the
> pseudopotential. Only *differences* (between structures, volumes, phases) are
> physical. We exploit that constantly from Chapter 3 on.

**A physical observable — the highest occupied level:**

```bash
grep "highest occupied" my.scf.out
```
```
     highest occupied level (ev):     6.2395
```
With 8 valence electrons filling 4 bands, this is the valence-band maximum.
(In Chapter 5 we turn this into a real band structure and read off the gap.)

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

Solutions: [`solutions/`](solutions/).

---

**Next:** [Chapter 1 — Plane-wave DFT & the SCF loop](../01-plane-wave-dft/),
where we open up the theory behind every line you just ran.
