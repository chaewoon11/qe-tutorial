---
id: nbse2-cdw
title: "MP1 — Charge density waves in monolayer NbSe₂"
sidebar_label: "MP1 · CDW in NbSe₂"
sidebar_position: 1
---

# Mini-project 1 — Charge density waves in monolayer NbSe₂

import Admonition from '@theme/Admonition';

This is the first **mini-project**: not a new technique, but a small piece of
*research* that chains together tools you already have — metals and smearing
([Ch. 3](../chapters/brillouin-zone.md)), DFPT phonons
([Ch. 10–11](../chapters/phonon-dispersion.md)), structural relaxation
([Ch. 4](../chapters/relaxation.md)), and charge densities
([Ch. 7](../chapters/charge.md)) — to reproduce a real many-body phenomenon: a
**charge density wave (CDW)** in monolayer 1H-NbSe₂.

The guiding idea is the one piece of new physics here:

> **A phonon dispersion can tell you that a crystal is unstable.** If a phonon
> frequency comes out **imaginary** (ω² < 0) at some wavevector **q**, the
> lattice *wants* to spontaneously distort with that periodicity. Following that
> soft mode down to a new, lower-energy structure **is** the CDW.

:::note Reference data
The phonon and supercell runs here are **HPC-grade** (a metal needs a dense
**k**-grid and a full DFPT **q**-grid). The inputs below are complete and runnable.
This page presents the **theory and the full workflow**; the quantitative results
(soft-mode frequency, condensation energy, charge-modulation figure) are being
finalized on the cluster and will be filled in as each stage completes. Spots
awaiting numbers are marked **[result to follow]**.
:::

---

## 1. The Peierls instability — why a 1D metal cannot stay uniform

Start with the cleanest case, a **1D chain** of atoms spaced by $a$ with one
electron per site (half-filled band). The Fermi wavevector is $k_F = \pi/2a$, so
the band is filled up to $\pm k_F$.

Now impose a small **periodic lattice distortion (PLD)** that doubles the period —
a displacement wave with wavevector $q = 2k_F = \pi/a$. Doubling the cell folds
the band at $\pm k_F$, and the electron–ion coupling opens a **gap** $2\Delta$
right there — exactly at the Fermi level. Every occupied state just below $k_F$ is
pushed **down** in energy; every empty state just above is pushed up. The
electrons therefore **lower** their total energy.

How the two energies scale with the distortion amplitude $u$ is the crux:

$$
\Delta E_\text{elastic} \sim +\tfrac12 K u^2 ,
\qquad
\Delta E_\text{electronic} \sim -\,u^2 \ln\!\frac{W}{u} ,
$$

where $K$ is the lattice stiffness and $W$ a bandwidth scale. The elastic cost is
a plain parabola, but the electronic gain carries a **logarithm** — and for small
$u$, $u^2\ln(1/u)$ always beats $u^2$. So the total energy is **lowered for any
nonzero $u$**: the undistorted 1D metal is *never* the ground state at $T=0$. This
is **Peierls' theorem**. The driving force is that the 1D Fermi "surface" (two
points $\pm k_F$) is **perfectly nested** by the single vector $2k_F$.

## 2. From 1D to a real crystal: nesting, susceptibility, and the Kohn anomaly

In 2D and 3D there is no automatic nesting, but the same machinery applies through
the **electronic susceptibility** (Lindhard response)

$$
\chi(\mathbf q) = \sum_{\mathbf k}
\frac{f(\varepsilon_{\mathbf k}) - f(\varepsilon_{\mathbf k+\mathbf q})}
     {\varepsilon_{\mathbf k+\mathbf q} - \varepsilon_{\mathbf k}} .
$$

If the Fermi surface has **parallel sheets** connected by one wavevector
$\mathbf q_\text{CDW}$ (a nesting vector), many $\mathbf k$-pairs contribute in
phase and $\chi(\mathbf q)$ develops a **peak** there. The electron–phonon
coupling then **renormalizes** the phonon frequency (Fröhlich / RPA):

$$
\omega^2_{\text{ren}}(\mathbf q)
   = \omega^2_{\text{bare}}(\mathbf q) - 2\,\omega_{\text{bare}}(\mathbf q)\,
     \frac{|g_{\mathbf q}|^2}{\hbar}\,\chi(\mathbf q) .
$$

A dip in $\omega(\mathbf q)$ at the nesting vector is a **Kohn anomaly**. Push the
coupling $|g_{\mathbf q}|^2\chi(\mathbf q)$ hard enough and the right-hand side
goes **negative**:

$$
\boxed{\;\omega^2_{\text{ren}}(\mathbf q_\text{CDW}) < 0
   \;\Longrightarrow\; \text{imaginary frequency} \;\Longrightarrow\;
   \text{lattice instability at } \mathbf q_\text{CDW}\;}
$$

The beauty of DFPT is that **we never build this model**. `ph.x` computes the
*fully renormalized* dynamical matrix directly, so a soft mode appears
automatically as an imaginary frequency. We just read it off the dispersion.

<Admonition type="tip" title="Why smearing (degauss) controls the CDW">
The susceptibility $\chi(\mathbf q)$ depends on how *sharp* the Fermi surface is.
The electronic smearing <code>degauss</code> acts like an electronic
**temperature**: a large <code>degauss</code> blurs the Fermi surface, weakens the
nesting peak, and the soft mode can disappear entirely — exactly as a real CDW
melts above its transition temperature $T_\text{CDW}$. So the depth of the soft
mode is genuinely <code>degauss</code>- and <strong>k</strong>-grid-dependent.
This is not a numerical nuisance; it <em>is</em> the physics, and we test it
explicitly below.
</Admonition>

## 3. The material: 1H-NbSe₂ and its 3×3 CDW

Bulk 2H-NbSe₂ is a classic CDW superconductor; a **single 1H layer** keeps the
CDW. The unit cell is one Nb in **trigonal-prismatic** coordination by six Se
(space group $P\bar 6m2$) — three atoms, a metal (Nb $4d^1$ gives a half-filled
$d_{z^2}$ band crossing $E_F$).

Experimentally the CDW is a **3×3** reconstruction. A $3\times3$ superstructure
folds wavevectors at multiples of $\tfrac13$ of the reciprocal lattice back to
$\Gamma$, so the instability we are hunting sits at

$$
\mathbf q_\text{CDW} = \left(\tfrac13, 0\right)
   \quad\text{— i.e. two-thirds of the way from } \Gamma \text{ to } M .
$$

That is why the phonon **q**-grid below is **6×6×1**: the point $(\tfrac13,0)$ is
$(\tfrac{2}{6},0)$, an *exact* grid point, so we sample the CDW wavevector without
interpolation.

The relaxed cell (Stage 0 below) comes out at **a₀ = 3.4736 Å**, with the Se
planes at **z = ±1.681 Å** — in line with PBE values for the monolayer.

---

## 4. Workflow

Four stages, each reusing the previous one's output — the professional discipline
of *converge → relax → reuse*:

| Stage | What | Tool | Cost |
|---|---|---|---|
| 0 | relax the undistorted 1H cell | `pw.x` vc-relax | light |
| 1 | phonon dispersion → find the soft mode | `pw.x`→`ph.x`→`q2r.x`→`matdyn.x` | **heavy** |
| 2 | freeze the soft mode into a 3×3 supercell, relax | `pw.x` relax | medium |
| 3 | visualize the charge modulation | `pp.x` (+ DOS) | light |

### Stage 0 — relax the high-symmetry cell

A metal, so we use Marzari–Vanderbilt cold smearing; `cell_dofree='2Dxy'` relaxes
the in-plane lattice constant while the vacuum spacing is held fixed, and
`assume_isolated='2D'` removes the spurious field between periodic images.

```fortran title="code/advanced/04-nbse2-cdw/inputs/nbse2.vc-relax.in"
&control
    calculation   = 'vc-relax'
    prefix        = 'nbse2'
    outdir        = './out'
    pseudo_dir    = '../../pseudos'
    verbosity     = 'high'
    forc_conv_thr = 1.0d-5
    etot_conv_thr = 1.0d-6
/
&system
    ibrav           = 0
    nat             = 3
    ntyp            = 2
    ecutwfc         = 60.0
    ecutrho         = 480.0
    assume_isolated = '2D'
    occupations     = 'smearing'
    smearing        = 'mv'
    degauss         = 0.02
/
&electrons
    conv_thr    = 1.0d-18
    mixing_beta = 0.7
/
&ions
    ion_dynamics = 'bfgs'
/
&cell
    cell_dynamics = 'bfgs'
    cell_dofree   = '2Dxy'
/
ATOMIC_SPECIES
  Nb 92.90638 Nb.pbe-spn-kjpaw_psl.1.0.0.UPF
  Se 78.971   Se.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
   3.480000000   0.000000000   0.000000000
  -1.740000000   3.013371920   0.000000000
   0.000000000   0.000000000  22.000000000
ATOMIC_POSITIONS crystal
  Nb  0.333333333  0.666666667   0.000000000
  Se  0.666666667  0.333333333   0.076136364
  Se  0.666666667  0.333333333  -0.076136364
K_POINTS automatic
  16 16 1 0 0 0
```

### Stage 1 — phonon dispersion of the undistorted cell

This is the heart of the project. We run a self-consistent calculation on the
relaxed cell with a **dense 18×18×1 k-grid** (a metal needs it), then a full DFPT
phonon calculation on a **6×6×1 q-grid**, interpolate the force constants with
`q2r.x`, and plot the dispersion with `matdyn.x` along Γ–M–K–Γ. NbSe₂ is a
**metal**, so there is *no* LO–TO splitting and *no* `epsil` — unlike the polar
insulators in the core chapters.

```fortran title="code/advanced/04-nbse2-cdw/inputs/nbse2.scf.in"
&control
    calculation = 'scf'
    prefix      = 'nbse2'
    outdir      = './out'
    pseudo_dir  = '../../pseudos'
    verbosity   = 'high'
/
&system
    ibrav           = 0
    nat             = 3
    ntyp            = 2
    ecutwfc         = 60.0
    ecutrho         = 480.0
    assume_isolated = '2D'
    occupations     = 'smearing'
    smearing        = 'mv'
    degauss         = 0.01
/
&electrons
    conv_thr    = 1.0d-16
    mixing_beta = 0.7
/
ATOMIC_SPECIES
  Nb 92.90638 Nb.pbe-spn-kjpaw_psl.1.0.0.UPF
  Se 78.971   Se.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
   3.473586449   0.000000000   0.000000000
  -1.736793225   3.008089368   0.000000000
   0.000000000   0.000000000  22.000000000
ATOMIC_POSITIONS crystal
  Nb  0.333333333  0.666666667   0.000000000
  Se  0.666666667  0.333333333   0.076425751
  Se  0.666666667  0.333333333  -0.076425751
K_POINTS automatic
  18 18 1 0 0 0
```

```fortran title="code/advanced/04-nbse2-cdw/inputs/nbse2.ph.in"
phonon dispersion of undistorted 1H-NbSe2 (metal, no LO-TO)
&inputph
  prefix   = 'nbse2'
  outdir   = './out'
  fildyn   = 'nbse2.dyn'
  fildvscf = 'dvscf'
  ldisp    = .true.
  nmix_ph  = 12
  tr2_ph   = 1.0d-16
  nq1 = 6, nq2 = 6, nq3 = 1
/
```

```fortran title="code/advanced/04-nbse2-cdw/inputs/nbse2.q2r.in"
&input
  fildyn = 'nbse2.dyn'
  flfrc  = 'nbse2.fc'
  zasr   = 'simple'
/
```

```fortran title="code/advanced/04-nbse2-cdw/inputs/nbse2.matdyn.in"
&input
  flfrc            = 'nbse2.fc'
  asr              = 'simple'
  flfrq            = 'nbse2.freq'
  q_in_band_form   = .true.
  q_in_cryst_coord = .true.
/
4
  0.000000  0.000000  0.0   80   ! Gamma
  0.500000  0.000000  0.0   50   ! M   (q_CDW = (1/3,0) lies 2/3 along here)
  0.333333  0.333333  0.0   80   ! K
  0.000000  0.000000  0.0    1   ! Gamma
```

**What to look for.** Along Γ→M one branch dives toward zero and crosses into
**imaginary** frequencies (plotted as negative) near $(\tfrac13,0)$ — the soft
mode. `matdyn.x` writes negative numbers there; that dip is the Kohn anomaly going
unstable, the DFT fingerprint of the CDW.

:::info Result to follow
Soft-mode frequency at $\mathbf q_\text{CDW}=(\tfrac13,0)$ (expected **imaginary**,
plotted negative) and the `degauss` sensitivity test will be reported here once the
DFPT run and the dispersion plot are complete.
:::

### Stage 2 — freeze the soft mode into a 3×3 supercell

Take the soft-mode **eigenvector** (the atomic displacement pattern) at
$\mathbf q_\text{CDW}$, build a commensurate **3×3 supercell** (27 atoms), displace
the atoms along that pattern, and **relax**. If the physics is right, the distorted
structure settles into a lower energy than the perfect 3×3 cell:

$$
\Delta E_\text{CDW} = E_\text{distorted} - E_\text{undistorted} < 0 ,
$$

the **condensation energy** per formula unit.

:::info Result to follow
The condensation energy $\Delta E_\text{CDW}$ (meV/formula unit) and the relaxed
distortion pattern — the characteristic Nb triangular clustering — will be reported
here once the 3×3 supercell relaxation completes.
:::

### Stage 3 — see the charge modulation

Finally, `pp.x` on the relaxed 3×3 cell gives the valence charge density; its
periodic modulation across the supercell **is** the "charge density wave." A DOS
comparison (undistorted vs distorted) shows the **pseudogap** the distortion opens
at $E_F$ — the energy-lowering mechanism from §1, now in a real 2D crystal.

:::info Result to follow
The charge-density figure across the 3×3 supercell and the DOS comparison
(undistorted vs distorted, showing the pseudogap at $E_F$) will be added here.
:::

---

## 5. Why this is the right capstone

Every ingredient is something you already met, now doing real work together:

- **smearing/metals** (Ch. 3) — without a metal there is no Fermi-surface nesting;
- **DFPT phonons** (Ch. 10–11) — the soft mode is the entire diagnostic;
- **relaxation** (Ch. 4) — the freeze-in needs a supercell BFGS relax;
- **charge density** (Ch. 7) — to *see* the wave.

And it connects directly to research you might do next: the **pre-ASR** force
constants computed here are exactly what feed electron–phonon (superconducting
$T_c$) and transport calculations, so getting the phonons well converged — not
just acoustic-sum-rule-patched — matters.

## Exercises

1. **Read the instability.** On the dispersion, identify the soft branch and the
   exact **q** where it is most negative. Does it match $(\tfrac13,0)$? What 2D
   superstructure does that wavevector imply?
2. **Melt the CDW with `degauss`.** Re-run Stage 1 with `degauss = 0.03` and
   `0.005`. Show that the soft mode shrinks (or vanishes) as `degauss` grows.
   Relate this to the experimental $T_\text{CDW}$.
3. **k-grid convergence.** The soft-mode depth depends on the **k**-grid feeding
   `ph.x`. Re-run with 12×12×1 and 24×24×1 scf grids. Which way does the mode move?
4. **Condensation energy.** From Stage 2, compute $\Delta E_\text{CDW}$ per formula
   unit. Is it of order a few meV? Compare to $k_B T_\text{CDW}$.
5. **No LO–TO here.** Explain why this calculation needs neither `epsil` nor
   `loto_2d`, even though the polar 2D crystal in [A3](../advanced/bilayer-hbn.md)
   did. (Hint: metal vs insulator.)

---

This closes the loop on the whole tutorial: the ground-state, electronic-structure,
and lattice-dynamics tools you built up, aimed together at a genuine condensed-matter
phenomenon.
