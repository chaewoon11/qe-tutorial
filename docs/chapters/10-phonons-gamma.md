---
title: "Chapter 10 — DFPT & Γ phonons"
sidebar_label: "10 · DFPT & Γ phonons"
sidebar_position: 11
---

# Chapter 10 — DFPT & Γ phonons

So far everything has been *static* — electrons relaxing in a fixed (or
optimizing) lattice. **Phonons** are the quantized vibrations of that lattice,
and the standard way to compute them in QE is **density-functional perturbation
theory (DFPT)** with `ph.x`. This chapter does the zone-centre (Γ) phonons of
GaAs, gets its Born effective charges and dielectric tensor (the ingredients for
LO–TO splitting), and runs the job on an HPC cluster.

**Prerequisites:** [Chapter 4](./04-relaxation.md) (forces) and a converged SCF.

---

## 1. From the harmonic approximation to the dynamical matrix

The nuclei sit near equilibrium positions $\mathbf{R}_I^0$ and vibrate about
them. Write the instantaneous position of nucleus $I$ as its equilibrium position
plus a small displacement,

$$
\mathbf{R}_I = \mathbf{R}_I^0 + \mathbf{u}_I,
$$

with Cartesian components $u_{I\alpha}$ ($\alpha = x,y,z$). Expand the
Born–Oppenheimer energy — the energy surface DFT gives us — as a Taylor series in
those displacements:

$$
E(\{\mathbf{u}\}) = E_0
+ \sum_{I\alpha}\left.\frac{\partial E}{\partial u_{I\alpha}}\right|_0 u_{I\alpha}
+ \frac12\sum_{I\alpha,J\beta} C_{I\alpha,J\beta}\,u_{I\alpha}u_{J\beta}
+ \mathcal{O}(u^3).
$$

At equilibrium the forces vanish, $\partial E/\partial u_{I\alpha}|_0 = -F_{I\alpha} = 0$,
so the **linear term drops out**. The **harmonic approximation** keeps only the
quadratic term, whose coefficients are the **interatomic force constants**

$$
C_{I\alpha,J\beta} = \left.\frac{\partial^2 E}{\partial u_{I\alpha}\,\partial u_{J\beta}}\right|_0 .
$$

The restoring force is then linear in the displacements (Hooke's law),
$F_{I\alpha} = -\sum_{J\beta} C_{I\alpha,J\beta}\,u_{J\beta}$, and Newton's
equations for the nuclei read

$$
M_I\,\ddot{u}_{I\alpha} = -\sum_{J\beta} C_{I\alpha,J\beta}\,u_{J\beta}.
$$

In a periodic crystal we seek **normal-mode (phonon) solutions** of Bloch form —
a wave of wavevector $\mathbf{q}$ and frequency $\omega$,

$$
u_{I\alpha}(t) = \frac{1}{\sqrt{M_I}}\,e_{I\alpha}(\mathbf{q})\,
e^{i(\mathbf{q}\cdot\mathbf{R}_I^0 - \omega t)} .
$$

Substituting into Newton's equations turns them into an **eigenvalue problem**,

$$
\omega^2(\mathbf{q})\,\mathbf{e}(\mathbf{q}) = D(\mathbf{q})\,\mathbf{e}(\mathbf{q}),
\qquad
D_{I\alpha,J\beta}(\mathbf{q}) = \frac{1}{\sqrt{M_I M_J}}
\sum_{\mathbf{R}} C_{I\alpha,J\beta}(\mathbf{R})\,e^{i\mathbf{q}\cdot\mathbf{R}},
$$

where $D(\mathbf{q})$ is the **dynamical matrix** (the mass-weighted Fourier
transform of the force constants). Diagonalizing this $3N_\text{at}\times3N_\text{at}$
Hermitian matrix gives the squared frequencies $\omega_{\mathbf{q}\nu}^2$
(eigenvalues) and the polarization vectors $\mathbf{e}_{\mathbf{q}\nu}$
(eigenvectors). A 2-atom cell like GaAs has $3\times2 = 6$ branches at each
$\mathbf{q}$: 3 acoustic + 3 optical.

So the *entire* harmonic phonon problem reduces to one ingredient — the force
constants $C = \partial^2 E/\partial u\,\partial u$. The rest of this chapter is
how DFPT computes them.

## 2. What the second derivative actually requires

Start from the force on an atom, which the Hellmann–Feynman theorem (Chapter 4)
gives the force as a ground-state expectation value, *without* differentiating
the wavefunctions:
$F_{I\alpha} = -\int n\,\partial_{I\alpha} v_\text{ext}\,d\mathbf{r} + F^\text{ion}_{I\alpha}$.
Differentiate **once more** and the density itself must now respond:

$$
\frac{\partial^2 E}{\partial u_{I\alpha}\partial u_{J\beta}} =
\int \frac{\partial n(\mathbf{r})}{\partial u_{J\beta}}\,
\frac{\partial v_\text{ext}(\mathbf{r})}{\partial u_{I\alpha}}\,d\mathbf{r}
+ \int n(\mathbf{r})\,
\frac{\partial^2 v_\text{ext}(\mathbf{r})}{\partial u_{I\alpha}\partial u_{J\beta}}\,d\mathbf{r}
+ \frac{\partial^2 E_\text{ion–ion}}{\partial u_{I\alpha}\partial u_{J\beta}} .
$$

Everything here is known from the ground state *except one object*: the
**first-order density response** $\partial n/\partial u_{J\beta}$ — how the
electrons rearrange when you nudge an atom. Computing that response *is* DFPT.

## 3. Linear response: the self-consistent Sternheimer equation

The density is built from the occupied Kohn–Sham orbitals, so its first-order
change comes from the first-order change of those orbitals,

$$
\Delta n(\mathbf{r}) = 4\,\mathrm{Re}\sum_{v}^{\text{occ}} \psi_v^*(\mathbf{r})\,\Delta\psi_v(\mathbf{r}),
$$

and each $\Delta\psi_v$ obeys a linear equation obtained by perturbing the
Kohn–Sham equation — the **Sternheimer equation**:

$$
(\hat{H}_\text{KS} - \varepsilon_v)\,\lvert\Delta\psi_v\rangle
= -\,\hat{P}_c\,\Delta v_\text{KS}\,\lvert\psi_v\rangle ,
$$

where $\hat{P}_c = 1-\sum_v|\psi_v\rangle\langle\psi_v|$ projects onto the empty
(conduction) subspace. The beauty: this needs **only the occupied orbitals** — no
sum over unoccupied states. The catch: the perturbing potential is itself
**self-consistent**, because moving the atom shifts the density, which shifts the
Hartree and xc potentials:

$$
\Delta v_\text{KS}(\mathbf{r}) = \underbrace{\Delta v_\text{ext}(\mathbf{r})}_{\text{the bare nudge}}
+ \int \frac{\Delta n(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}\,d\mathbf{r}'
+ \left.\frac{dv_\text{xc}}{dn}\right|_{n(\mathbf{r})}\!\!\Delta n(\mathbf{r}) .
$$

So DFPT runs a **second self-consistent loop** — exactly like the ground-state
SCF of Chapter 1, but now for the *response* $\Delta n$ instead of $n$: guess
$\Delta n$ → build $\Delta v_\text{KS}$ → solve the Sternheimer equation for
$\Delta\psi_v$ → rebuild $\Delta n$ → repeat. **`tr2_ph` is the `conv_thr` of
this response loop** (§ note below).

## 4. Why DFPT — and any **q** for free

Two payoffs make this worth it over brute-force finite differences ("frozen
phonons", where you displace atoms in a supercell and re-run SCF):

- **The 2n+1 theorem.** Knowing the wavefunctions to order $n$ determines the
  energy to order $2n+1$. The *first*-order response $\Delta\psi$ ($n=1$) thus
  fixes the energy to *third* order — enough for the second derivative (and even
  third-order anharmonic constants). You never need the second-order
  wavefunctions.
- **Each q is independent.** A displacement wave $u_{I\alpha}\propto e^{i\mathbf{q}\cdot\mathbf{R}_I}$
  is *monochromatic*: it couples a Bloch state at $\mathbf{k}$ only to $\mathbf{k}+\mathbf{q}$.
  So the response at wavevector $\mathbf{q}$ can be computed **in the primitive
  cell**, one q at a time, at **arbitrary q** — no supercell. (Frozen phonons
  need a supercell commensurate with $\mathbf{q}$, which explodes for general q.)

`ph.x` implements all of this. The canonical reference is Baroni, de Gironcoli,
Dal Corso & Giannozzi, *Rev. Mod. Phys.* **73**, 515 (2001).

:::note Why `tr2_ph = 1e-16` ≈ an SCF `conv_thr` of `1e-8`
`ph.x` has its own convergence threshold, **`tr2_ph`**, and it looks far tighter
than a normal SCF's `conv_thr` — but they measure different things. DFPT works
with a **perturbation**: the small quantity is the *first-order* change of the
wavefunctions/density, $\Delta^{(1)}$. `tr2_ph` is a threshold on a quantity that
scales like its **square**, $|\Delta^{(1)}|^2$. So `tr2_ph = 10^{-16}` controls
the first-order response to $\sqrt{10^{-16}} = 10^{-8}$ — the same *linear*
accuracy a ground-state SCF gets from `conv_thr = 10^{-8}`. That squared
convention is why phonon thresholds carry "extra" zeros; `1e-16` here is a normal,
well-converged setting, not an extreme one.
:::

## 5. Hands-on: Γ phonons of GaAs

The recipe is **converged SCF → `ph.x`**. The phonon run reads the SCF
wavefunctions from `outdir`, so both steps share the same directory. Because GaAs
is **polar**, we set `epsil = .true.` to also compute the **dielectric tensor**
and **Born effective charges** — the ingredients for the LO–TO splitting
(Chapter 11).

```fortran title="code/10-phonons-gamma/gaas.ph.gamma.in"
GaAs Gamma-point phonons (LO-TO via Born effective charges)
&inputph
    prefix   = 'gaas'
    outdir   = './out'
    fildyn   = 'gaas.dyn'
    tr2_ph   = 1.0d-16
    nmix_ph  = 12         ! mix more response densities -> faster convergence
    epsil    = .true.     ! dielectric tensor + Born charges (polar material)
    trans    = .true.     ! compute the phonons (dynamical matrix)
/
0.0 0.0 0.0               ! the q-point: Γ
```

The SCF input is the usual GaAs insulator (`occupations='fixed'`), but **phonons
demand tighter convergence than the ground state did**. The pseudopotentials only
*suggest* `ecutwfc = 33` Ry / `ecutrho = 134` Ry (the larger of Ga's 33/134 and
As's 20/103) — a floor for *energies*. Phonons are second energy derivatives and
are far more sensitive, so we use **`ecutwfc = 70`, `ecutrho = 560`, a 12×12×12
k-grid, and `conv_thr = 1e-20`** (we'll justify those numbers in §3 with the
acoustic modes as the convergence gauge). This is the first **genuinely
expensive** step in the tutorial, so we ran it on **Nurion's `debug` queue**
(KNL) rather than a workstation:

```bash
# in the PBS job (debug queue, 1 KNL node):
mpirun -np 68 pw.x -npool 4 < gaas.scf.in      > gaas.scf.out
mpirun -np 68 ph.x -npool 4 < gaas.ph.gamma.in > gaas.ph.out
```

It finished in well under ten minutes on one node. `nmix_ph = 12` (mixing more
response-density histories, like a larger `mixing_ndim`) noticeably cuts the
number of DFPT iterations.

:::tip Run phonons on HPC
DFPT scales steeply with cell size, q-points, and the number of irreducible
perturbations. Even this 2-atom Γ job is minutes on 68 cores; a full dispersion
(Chapter 11) is many such jobs. Use a cluster — see [`docs/`](https://github.com/chaewoon11/qe-tutorial/tree/master/reference/hpc-guides)
for our Nurion / Stampede3 / Perlmutter notes.
:::

## 6. Reading the results

**Frequencies at Γ** (`grep freq gaas.ph.out`):

```text
     freq (  1 -  3) =        -3.5  [cm-1]   --> T_2  G_15 P_4   I+R
     freq (  4 -  6) =       269.7  [cm-1]   --> T_2  G_15 P_4   I+R
```

- Modes 1–3 are the **acoustic** branches. They must go to **zero** at Γ (a rigid
  translation of the whole crystal costs no energy). Here they sit at **−3.5 cm⁻¹**
  (the sign just means a slightly negative eigenvalue of the raw matrix) — small,
  as a converged calculation should be. More on that below.
- Modes 4–6 are the triply-degenerate **optical** mode, the **TO** frequency
  **269.7 cm⁻¹** (8.10 THz) — right on the experimental ~268 cm⁻¹. The `I+R`
  label says it is both infrared- and Raman-active.

:::warning The acoustic modes are your convergence gauge
`ph.x` prints the **raw** dynamical-matrix frequencies — it does *not* apply the
acoustic sum rule. So the Γ acoustic frequencies, which translational invariance
says must be exactly **0**, are a free diagnostic of how converged the *whole*
calculation is. They depend on the **k-grid and cutoffs**, not on the SCF
`conv_thr` — watch what happens as we tighten them:

| `ecutwfc` / `ecutrho` | k-grid | acoustic at Γ (raw) |
|---|---|---|
| 50 / 400 | 8×8×8 | 13.7 cm⁻¹ — **under-converged** |
| 50 / 400 | 8×8×8, `conv_thr` 1e-12 → 1e-20 | 13.4 cm⁻¹ — *barely moved* |
| **70 / 560** | **12×12×12** | **3.5 cm⁻¹ — converged** |

The pseudopotentials' "suggested" 33/134 Ry, and even our ground-state 50/400, are
*too loose for phonons*: second derivatives need a denser k-grid and higher
cutoffs. Note the middle row — tightening the SCF `conv_thr` from 1e-12 to 1e-20
did essentially nothing, which tells you the residual was a **basis/k-sampling**
error, not an SCF one.

And the key point: **the acoustic sum rule (ASR) merely forces those three modes
to exactly 0** — a symmetry constraint imposed by hand, not a cure. A large raw
violation (the 13.7 cm⁻¹ first row) means the *rest* of the spectrum — the
optical mode and every finite-**q** frequency in Chapter 11 — is under-converged
too, by an amount you can no longer see once ASR has zeroed the acoustic branch.
So converge first (here to 3.5 cm⁻¹), *then* let ASR clean up the residual.
:::

:::note Why this matters downstream
This is not just cosmetic. Whenever phonons feed a **later** calculation —
electron–phonon coupling in EPW, polaron formation energies, superconducting
$T_c$ via the Eliashberg function $\alpha^2F(\omega)$ and $\lambda$, transport,
or thermal conductivity — what those codes consume is the **raw dynamical
matrices / interatomic force constants**, i.e. the *pre-ASR* phonon data across
the whole Brillouin zone. ASR only repairs the three translational modes at Γ; it
cannot fix force constants that are under-converged everywhere else. An
under-converged phonon spectrum propagates straight into $\lambda$, $T_c$, polaron
binding energies, and the rest. So for any production electron–phonon workflow,
converge the phonons *properly first* (cutoffs, k- and q-grids, thresholds) — a
small pre-ASR acoustic frequency is your cheapest sign that you have.
:::

**Dielectric tensor** (`epsil`):

```text
     Dielectric constant in cartesian axis
        17.837   0.000   0.000
         0.000  17.837   0.000
         0.000   0.000  17.837
```

Isotropic (scalar) as required by cubic symmetry. $\varepsilon_\infty = 17.8$ is
still above the experimental **10.9** — a residue of the **PBE band-gap problem**
from Chapter 5 ($\varepsilon_\infty$ grows as the gap shrinks, so PBE's too-small
gap over-polarizes the crystal). Tellingly, the loose 50/400 run gave a wildly
inflated 29.3; converging the basis nearly halved it, another reminder that the
under-converged numbers were not to be trusted.

**Born effective charges** $Z^*$ (how much polarization an atom's displacement
creates):

```text
   atom 1 Ga   Mean Z*:  +2.160     (with acoustic sum rule applied)
   atom 2 As   Mean Z*:  -2.160
```

$Z^* = \pm2.16$ matches the experimental value almost exactly — the fingerprint
of a **polar** crystal, and what splits the LO mode from the TO mode. Note the
**acoustic sum rule** at work: the raw $Z^*$ (Ga $+2.140$, As $-2.181$) sum to
$-0.041$ instead of zero (a small residual, vs. the loose run's $-0.247$); QE's
`asr` symmetrizes them to $\pm2.160$ so a rigid translation produces no net force.

The TO frequency, $Z^*$, and $\varepsilon_\infty$ are exactly the three
ingredients needed to build the **LO–TO splitting** and the full phonon
dispersion — which is Chapter 11.

## Exercises

1. **Acoustic sum rule.** In the loose run the acoustic modes came out at
   13.7 cm⁻¹ instead of 0 and the raw $Z^*$ summed to $-0.247$; converging the
   basis/k-grid brought these to −3.5 cm⁻¹ and $-0.041$. What single physical
   principle do both violations break, and why is converging the calculation the
   real fix rather than just applying ASR?
2. **Threshold intuition.** If `tr2_ph = 1e-16` corresponds to a linear accuracy
   of $10^{-8}$, what `tr2_ph` would you choose to match an SCF converged to
   `conv_thr = 1e-10`? Why is going much tighter usually wasted effort?
3. **The dielectric red flag.** GaAs's computed $\varepsilon_\infty = 17.8$ vs.
   experiment 10.9 (and 29.3 in the under-converged run!). Connect this to the
   Chapter 5 band gap. Would a hybrid functional raise or lower
   $\varepsilon_\infty$?
4. **Mode counting.** Why does GaAs have exactly 6 phonon modes at Γ, and why are
   the 3 optical modes degenerate? What would change for a 4-atom cell?

<details>
<summary>Show solutions</summary>

**1. Acoustic sum rule.** Both break **translational invariance**: rigidly
shifting the entire crystal must cost zero energy (→ acoustic frequencies = 0 at
Γ) and produce zero net force (→ $\sum_I Z^*_I = 0$). Numerical incompleteness
(finite grids/cutoffs) violates it; imposing the **acoustic sum rule** restores
exact translational invariance, zeroing the acoustic modes and making the $Z^*$
sum vanish. But ASR only *enforces the symmetry* — it does not improve the
underlying convergence. The loose run's 13.7 cm⁻¹ violation is large; ASR would
zero it but leave the optical and finite-**q** frequencies under-converged. The
honest fix is tighter cutoffs/k-grid (which got us to 3.5 cm⁻¹), with ASR as the
final touch-up — not a substitute for converging.

**2. Threshold intuition.** Since `tr2_ph` $\sim$ (linear accuracy)², matching
`conv_thr = 1e-10` means `tr2_ph = (1e-10)² = 1e-20`. Going much tighter than the
accuracy you actually need just burns DFPT iterations without changing the
frequencies — `1e-16`–`1e-18` is plenty for most work.

**3. The dielectric red flag.** $\varepsilon_\infty$ is dominated by virtual
valence→conduction transitions, whose strength scales roughly as $1/E_g^2$. PBE's
far-too-small gap (0.5 eV vs. 1.42 eV) inflates $\varepsilon_\infty$. A hybrid
functional opens the gap toward experiment and would **lower**
$\varepsilon_\infty$ closer to 10.9.

**4. Mode counting.** Modes $= 3 \times N_\text{atoms} = 3\times2 = 6$. The three
optical modes are degenerate because the zincblende point group makes the three
Cartesian directions equivalent (a $T_2$ representation). A 4-atom cell would have
12 modes, and lower symmetry would generally lift some degeneracies.

</details>

---

**Next:** [Chapter 11 — Phonon dispersion](./11-phonon-dispersion.md): a q-grid
with `ph.x`, Fourier interpolation with `q2r.x` / `matdyn.x`, the acoustic sum
rule, and the LO–TO splitting that the Born charges above produce.
