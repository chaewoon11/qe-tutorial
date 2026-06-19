---
title: "Chapter 4 — Structural optimization"
sidebar_label: "4 · Structural optimization"
sidebar_position: 5
---

# Chapter 4 — Structural optimization

So far we fixed the geometry by hand (the experimental GaAs lattice constant).
But the whole point of DFT is often to *predict* a structure: relax the atoms and
the cell until they sit at the energy minimum. This chapter covers the forces
and stresses that drive that relaxation, the `relax` / `vc-relax` calculations
that perform it, and a classic result — PBE overestimates lattice constants.

**Prerequisites:** [Chapter 3](./03-brillouin-zone.md) (you should converge the
k-grid before trusting a relaxation).

---

## 1. Forces and stress: the gradients of the energy

A structure is in equilibrium when the total energy is stationary — no net force
on any atom and no residual stress on the cell. Both are *derivatives* of the
Kohn–Sham energy, and DFT gives them almost for free thanks to the
**Hellmann–Feynman theorem**: because the energy is variational, the force on an
atom is just the expectation value of the gradient of the potential, with no need
to differentiate the wavefunctions,

$$
\mathbf{F}_I = -\frac{\partial E}{\partial \mathbf{R}_I}
= -\Big\langle \psi \Big| \frac{\partial \hat{V}_\text{ext}}{\partial \mathbf{R}_I} \Big| \psi \Big\rangle .
$$

The analogous derivative with respect to a uniform strain $\varepsilon_{\alpha\beta}$
of the cell is the **stress tensor**

$$
\sigma_{\alpha\beta} = \frac{1}{\Omega}\frac{\partial E}{\partial \varepsilon_{\alpha\beta}} .
$$

To get the stress (needed for `vc-relax`) you must tell `pw.x` to compute it:

```fortran
&control
    calculation = 'vc-relax'
    tprnfor     = .true.    ! print forces
    tstress     = .true.    ! print the stress tensor
/
```

QE prints forces in Ry/bohr (`Forces acting on atoms`) and stress in both
kbar and Ry/bohr³ (`total   stress`). At the minimum both should vanish to within
your thresholds.

## 2. `relax` vs `vc-relax`

| Calculation | What moves | What it finds |
|---|---|---|
| `relax` | atomic positions only, **cell fixed** | internal coordinates at fixed cell shape/size |
| `vc-relax` | atomic positions **and** the cell | the full equilibrium structure (cell + atoms) |

You need the `&ions` namelist for `relax`, plus `&cell` for `vc-relax`:

```fortran
&ions
    ion_dynamics = 'bfgs'      ! quasi-Newton optimizer
/
&cell
    cell_dofree  = 'all'       ! which cell DOF may change (see below)
    press        = 0.0d0       ! target pressure (kbar)
/
```

The optimizer is **BFGS**, a quasi-Newton method: it uses the forces (and an
accumulated approximate Hessian) to step downhill toward the minimum, re-running
the SCF at each new geometry. `cell_dofree` restricts which cell degrees of
freedom may relax — `'all'` (everything), `'ibrav'` (keep the Bravais symmetry),
`'2Dxy'` (relax in-plane only, for slabs), `'volume'`, etc.

:::note GaAs has no internal coordinates
In zincblende the two atoms sit at $(0,0,0)$ and $(\tfrac14,\tfrac14,\tfrac14)$,
fixed by symmetry — so a plain `relax` would have nothing to do. Only the
**cell** (one number, the lattice constant) can relax. That makes GaAs a clean
`vc-relax` example; we'll demonstrate `relax`-style forces by displacing an atom
in the exercises.
:::

## 3. Convergence thresholds and a gotcha

Relaxation stops when the changes fall below thresholds you set:

| Threshold | Meaning | Typical |
|---|---|---|
| `forc_conv_thr` | max force on any atom | `1.0d-4` – `1.0d-3` Ry/bohr (≈ 1–25 meV/Å) |
| `etot_conv_thr` | energy change between steps | `1.0d-4` – `1.0d-5` Ry |
| (stress) | residual stress for `vc-relax` | controlled by the same force/energy criteria + `press_conv_thr` |

:::warning Re-converge after a `vc-relax`
A plane-wave basis is defined by `ecutwfc` at a *fixed* cell. When `vc-relax`
changes the cell volume, the effective basis changes, which introduces a
spurious **Pulay stress**. The fix is to run a **fresh `scf` (or a second
`vc-relax`) on the relaxed cell** to get a clean, consistent final energy and
stress. Always do this final single-point before quoting numbers.
:::

## 4. Hands-on: the equation of state of GaAs

The most transparent way to *see* the equilibrium structure is to scan it: run
`scf` at a series of lattice constants and plot the total energy. Each point is a
fixed-geometry calculation like Chapter 0, just with a different `celldm(1)`:

```bash
cd code/04-relaxation
for a in 10.20 10.40 10.60 10.70 10.80 10.90 11.00 11.20; do
  sed "s/celldm(1)     = .*/celldm(1)     = $a/" gaas.scf.in > /tmp/e.in
  echo "$a $(grep '^!' <(pw.x < /tmp/e.in) | tail -1 | awk '{print $5}')"
done
```

[`notebooks/04-equation-of-state.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/04-equation-of-state.ipynb)
fits these points to the **Birch–Murnaghan** equation of state, which gives the
equilibrium lattice constant $a_0$ and the bulk modulus $B_0$:

![GaAs equation of state](/img/04-eos.png)

The fit yields:

- **PBE equilibrium $a_0 \approx 5.756$ Å**, versus the experimental
  **5.653 Å** — PBE overestimates by **+1.8 %**, the systematic GGA error for
  semiconductor lattice constants.
- **Bulk modulus $B_0 \approx 60$ GPa** (experiment ≈ 75 GPa) — PBE
  *under*estimates the stiffness, the flip side of the too-large cell.

## 5. Hands-on: let `vc-relax` find it automatically

Scanning is instructive but wasteful — `vc-relax` walks straight to the minimum
by following the stress. The input
[`code/04-relaxation/gaas.vc-relax.in`](https://github.com/chaewoon11/qe-tutorial/blob/master/code/04-relaxation/gaas.vc-relax.in)
starts from the experimental lattice constant and relaxes the cell:

```bash
cd code/04-relaxation
pw.x < gaas.vc-relax.in > vc-relax.out
grep -A6 "Begin final coordinates" vc-relax.out
```

The optimizer drives the stress to zero and reports the relaxed cell. Here it
converges in just **3 BFGS steps** to a final cell volume of 320.4 bohr³, i.e. a
relaxed lattice constant of **5.748 Å** — the same ~5.75 Å minimum the
equation-of-state scan found, reached in a handful of steps instead of a full
sweep. (Remember §3: a final `scf` on this cell removes Pulay
stress before you quote an energy.)

## Exercises

1. **Read off the bulk modulus units.** The notebook reports $B_0$ in GPa from a
   fit in Ry/bohr³. Verify the conversion $1\ \mathrm{Ry/bohr^3} = 14710.5$ GPa
   makes the GaAs $B_0$ come out near the experimental ~75 GPa.
2. **Hellmann–Feynman in action.** Displace the As atom from
   $(\tfrac14,\tfrac14,\tfrac14)$ to $(0.26,0.25,0.25)$ and run a `relax`
   (`calculation='relax'`, add `&ions`). Where does it return to, and what does
   the initial force on As tell you?
3. **PBE vs experiment.** Your $a_0$ is ~1 % too large. Why does that matter for,
   say, a predicted phonon frequency or band gap? (Think about what later
   chapters compute *at* this lattice constant.)
4. **Why re-converge?** Run a single `scf` at the relaxed lattice constant and
   compare its total energy / pressure to the last step of the `vc-relax`. How
   big is the Pulay-stress artefact at `ecutwfc = 50`?

<details>
<summary>Show solutions</summary>

**1. Bulk modulus units.** With $B_0$ from the fit in Ry/bohr³, multiply by
14710.5 to get GPa. For GaAs the fit gives ≈ 60 GPa, in the right ballpark of the
experimental ~75 GPa (PBE *under*estimates $B_0$, hand in hand with the
over-large lattice constant).

**2. Hellmann–Feynman.** The displaced As feels a restoring force pointing back
toward $(\tfrac14,\tfrac14,\tfrac14)$; `relax` returns it there (to within
`forc_conv_thr`), because that symmetric site is the energy minimum and the force
is $-\partial E/\partial\mathbf{R}$. The magnitude of the *initial* force is the
slope of the energy along that displacement.

**3. PBE vs experiment.** Everything downstream is computed *at* the lattice
constant you choose. A 1 % too-large cell softens bonds: phonon frequencies come
out too low, and the band gap shifts (a larger lattice generally narrows the gap).
This is why some workflows use the *experimental* lattice constant, or a better
functional, when the property is lattice-sensitive.

**4. Why re-converge.** The single-point `scf` energy/pressure on the relaxed
cell differ slightly from the `vc-relax` final step because the basis was
defined at the *old* cell during relaxation (Pulay stress). At `ecutwfc = 50`
the residual stress is small but nonzero; pushing `ecutwfc` higher shrinks it.

</details>

---

**Next:** Chapter 5 — Band structure *(coming soon)*: from the SCF density to
`nscf` eigenvalues along a k-path, and reading off the (in)direct gap.
