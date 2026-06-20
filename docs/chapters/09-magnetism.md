---
title: "Chapter 9 — Spin & magnetism"
sidebar_label: "9 · Spin & magnetism"
sidebar_position: 10
---

# Chapter 9 — Spin & magnetism

Chapter 8 already turned on spin to get NiO's antiferromagnetism. This chapter
treats magnetism properly: how collinear spin-polarized DFT works and exactly how
QE reports a magnetic moment, the ferromagnet **bcc iron** as the worked example,
and what changes for **noncollinear** magnetism and **spin–orbit coupling**.

**Prerequisites:** [Chapter 8](./08-functionals.md) (spin was introduced there);
[Chapter 6](./06-dos.md) for the DOS.

---

## 1. Collinear spin-polarized DFT

With `nspin = 2`, the electron density splits into two spin channels, each built
from its own Kohn–Sham orbitals,

$$
n_\sigma(\mathbf{r}) = \sum_{n\mathbf{k}} f_{n\mathbf{k}\sigma}\,
\lvert \psi_{n\mathbf{k}\sigma}(\mathbf{r})\rvert^2, \qquad \sigma\in\{\uparrow,\downarrow\},
$$

and the exchange–correlation potential becomes **spin-dependent**,
$v_{\text{xc}}^{\sigma}(\mathbf{r}) = \delta E_\text{xc}[n_\uparrow,n_\downarrow]/\delta n_\sigma$,
so up and down electrons feel different potentials and can occupy differently —
the origin of spontaneous magnetism. The two natural fields are the total density
and the **magnetization density**:

$$
n(\mathbf{r}) = n_\uparrow + n_\downarrow, \qquad m(\mathbf{r}) = n_\uparrow - n_\downarrow .
$$

:::note How QE reports the magnetic moment
QE prints two cell integrals of $m(\mathbf{r})$, both in $\mu_B$/cell:

$$
M_\text{tot} = \int m(\mathbf{r})\,d\mathbf{r} = N_\uparrow - N_\downarrow,
\qquad
M_\text{abs} = \int \lvert m(\mathbf{r})\rvert\,d\mathbf{r}.
$$

**`total magnetization`** is the signed net moment ($N_\uparrow-N_\downarrow$
unpaired electrons, since $g\!\approx\!2$). **`absolute magnetization`**
integrates the *magnitude*, so it stays positive even when local moments cancel.
A ferromagnet has $M_\text{tot}=M_\text{abs}$; an antiferromagnet has
$M_\text{tot}=0$ but $M_\text{abs}>0$ (NiO in Chapter 8: 0 vs 3.08). The
**per-site** moment is the same integral over a sphere around each atom,
$m_I = \int_{|\mathbf{r}-\mathbf{R}_I|<R_I} m(\mathbf{r})\,d\mathbf{r}$.
:::

A metal turns ferromagnetic when the exchange energy gained by splitting the
spins outweighs the kinetic cost of moving electrons to higher states — the
**Stoner criterion**, $I\,g(E_F) > 1$, where $I$ is the exchange integral and
$g(E_F)$ the (nonmagnetic) density of states at the Fermi level. Iron, with its
narrow, high-DOS 3d bands, satisfies it comfortably.

## 2. Hands-on: bcc iron, a ferromagnet

Body-centred-cubic Fe is *the* itinerant ferromagnet. The input is a one-atom
metal with `nspin = 2` and a nonzero `starting_magnetization` to break spin
symmetry:

```fortran title="code/09-magnetism/fe.fm.scf.in"
&control
    calculation = 'scf'
    prefix      = 'fe'
    outdir      = './out'
    pseudo_dir  = '../pseudos'
    verbosity   = 'high'
/
&system
    ibrav = 3                       ! bcc
    celldm(1) = 5.4235              ! a = 2.87 Å
    nat = 1
    ntyp = 1
    ecutwfc = 80.0
    ecutrho = 640.0
    nbnd = 16
    nspin = 2                       ! collinear spin polarization
    starting_magnetization(1) = 0.7 ! initial guess (fraction of max)
    occupations = 'smearing'        ! metal
    smearing = 'mv'
    degauss = 0.02
/
&electrons
    mixing_beta = 0.3
    conv_thr = 1.0d-9
/
ATOMIC_SPECIES
  Fe 55.845 Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Fe 0.0 0.0 0.0
K_POINTS automatic
  12 12 12 0 0 0
```

Run it, and run a nonmagnetic (`nspin = 1`) version for comparison — the
nonmagnetic input is the same cell with `nspin = 1` and no
`starting_magnetization`:

```fortran title="code/09-magnetism/fe.nm.scf.in"
&control
    calculation = 'scf'
    prefix      = 'fe_nm'
    outdir      = './out_nm'
    pseudo_dir  = '../pseudos'
    verbosity   = 'high'
/
&system
    ibrav = 3
    celldm(1) = 5.4235
    nat = 1
    ntyp = 1
    ecutwfc = 80.0
    ecutrho = 640.0
    nbnd = 16
    nspin = 1
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.02
/
&electrons
    mixing_beta = 0.3
    conv_thr = 1.0d-9
/
ATOMIC_SPECIES
  Fe 55.845 Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Fe 0.0 0.0 0.0
K_POINTS automatic
  12 12 12 0 0 0
```

```bash
cd code/09-magnetism
pw.x < fe.fm.scf.in > fe.fm.out     # ferromagnetic
pw.x < fe.nm.scf.in > fe.nm.out     # nonmagnetic (nspin=1)
grep -E "total magnetization|^!" fe.fm.out
```

Results:

| | nonmagnetic | ferromagnetic | experiment |
|---|---|---|---|
| moment per Fe | 0 (imposed) | 2.23 μB | 2.22 μB |
| total energy | reference | −0.57 eV/atom | — |

The moment of **2.23 μB** lands right on the textbook iron value (the per-site
sphere integral gives 2.28 μB), and the ferromagnetic state is **lower in energy**
than the nonmagnetic one by **0.57 eV/atom** — the magnetic stabilization energy
that makes Fe a magnet at all. (Note `total magnetization` 2.23 vs `absolute
magnetization` 2.36: they nearly coincide because Fe is a ferromagnet, with a
little extra from the spin density between the atomic spheres.)

## 3. Spin-resolved density of states

The mechanism is clearest in the DOS. With `nspin = 2`, `dos.x` writes separate
**majority** (↑) and **minority** (↓) densities. A dense, spin-polarized `nscf`
(tetrahedra) feeds `dos.x`:

```fortran title="code/09-magnetism/fe.nscf.in"
&control
    calculation = 'nscf'
    prefix      = 'fe'
    outdir      = './out'
    pseudo_dir  = '../pseudos'
    verbosity   = 'high'
/
&system
    ibrav = 3
    celldm(1) = 5.4235
    nat = 1
    ntyp = 1
    ecutwfc = 80.0
    ecutrho = 640.0
    nbnd = 16
    nspin = 2
    starting_magnetization(1) = 0.7
    occupations = 'tetrahedra'
/
&electrons
    conv_thr = 1.0d-9
/
ATOMIC_SPECIES
  Fe 55.845 Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Fe 0.0 0.0 0.0
K_POINTS automatic
  16 16 16 0 0 0
```

```fortran title="code/09-magnetism/fe.dos.in"
&dos
    prefix = 'fe'
    outdir = './out'
    fildos = 'fe.dos'
    Emin = 8.0
    Emax = 25.0
    DeltaE = 0.02
/
```

Plotting the two channels back-to-back:

![bcc Fe spin-resolved DOS](/img/09-fe-dos.png)

The two spin channels are **rigidly shifted** relative to each other by the
**exchange splitting** Δ_ex: the majority 3d band sits lower and is almost fully
occupied, the minority band is pushed up and partly empty. The moment is just the
difference in filled states, $M = N_\uparrow - N_\downarrow$ — the area imbalance
between the two halves of the plot. This is the Stoner picture made visible.

## 4. Noncollinear magnetism & spin–orbit coupling

Collinear DFT forces every moment onto a single global axis (↑ or ↓). Real
materials can have moments pointing in **different directions** (canted spins,
spirals, skyrmions) and feel **spin–orbit coupling** (SOC), which ties the spin
direction to the lattice. For that you switch to a **noncollinear** description:
the orbitals become two-component **spinors** $\psi = (\psi^\uparrow, \psi^\downarrow)^T$,
and the magnetization is a **vector field**

$$
\mathbf{m}(\mathbf{r}) = \sum_{n\mathbf{k}} f_{n\mathbf{k}}\,
\psi_{n\mathbf{k}}^\dagger(\mathbf{r})\,\boldsymbol{\sigma}\,\psi_{n\mathbf{k}}(\mathbf{r}),
$$

with $\boldsymbol{\sigma}$ the Pauli matrices. In QE:

```fortran
&system
    noncolin   = .true.     ! vector magnetization (spinor wavefunctions)
    lspinorb   = .true.     ! add spin–orbit coupling (needs fully-relativistic pseudos)
    starting_magnetization(1) = 0.7
    angle1(1) = 90, angle2(1) = 0   ! initial moment direction (polar/azimuthal, degrees)
/
```

Two caveats make this much heavier than collinear: it **doubles the basis**
(spinors) and **requires fully-relativistic (`rel='full'`) pseudopotentials** for
SOC. What it buys you: **magnetocrystalline anisotropy** (the tiny energy
difference between moment along $[001]$ vs $[111]$ — why a magnet has easy axes),
SOC band splittings, and genuinely noncollinear spin textures. These are
specialist calculations; the collinear `nspin=2` above covers the large majority
of magnetism problems.

## Exercises

1. **Moment from the DOS.** From the spin-resolved DOS, the moment is the
   imbalance $N_\uparrow-N_\downarrow$ in occupied states. Confirm it equals the
   `total magnetization` QE printed.
2. **Kill the magnetism.** Set `starting_magnetization = 0` (with `nspin=2`).
   Does Fe stay magnetic or collapse to the nonmagnetic solution? What does that
   tell you about the role of the *initial* guess vs. the self-consistent result?
3. **Smearing and metals.** Why must bcc Fe use `occupations='smearing'`
   (Chapter 3), and how could too large a `degauss` corrupt the computed moment?
4. **Why noncollinear for anisotropy?** Explain why magnetocrystalline anisotropy
   energy *cannot* come out of a collinear `nspin=2` calculation, no matter how
   well converged.

<details>
<summary>Show solutions</summary>

**1. Moment from the DOS.** Integrating the majority and minority DOS up to $E_F$
gives $N_\uparrow$ and $N_\downarrow$; their difference equals the printed
`total magnetization` (≈ 2.23 μB). The moment *is* the spin-channel population
imbalance.

**2. Kill the magnetism.** With exactly zero starting magnetization the up and
down channels are identical by symmetry, and the SCF often stays stuck in the
(metastable) nonmagnetic solution — you must *break the symmetry* with a nonzero
initial guess to let it find the lower-energy ferromagnetic state. The converged
result is physical; the initial guess only decides which basin you fall into.

**3. Smearing and metals.** Fe is a metal with a Fermi surface, so a finite
k-grid needs smearing to integrate the discontinuous occupations (Chapter 3). Too
large a `degauss` over-smears the sharp 3d DOS at $E_F$, blurring the exchange
splitting and **underestimating the moment** — converge `degauss` together with
the k-grid.

**4. Why noncollinear for anisotropy.** Magnetocrystalline anisotropy comes
entirely from **spin–orbit coupling**, which couples the spin direction to the
lattice. Collinear `nspin=2` has no SOC and no notion of moment *direction*
relative to the crystal axes — every orientation is identical in energy. You need
`noncolin=.true., lspinorb=.true.` (and fully-relativistic pseudos) to get a
nonzero anisotropy.

</details>

---

**Next:** [Chapter 10 — DFPT & Γ phonons](./10-phonons-gamma.md): lattice
dynamics from density-functional perturbation theory, starting at the zone centre.
