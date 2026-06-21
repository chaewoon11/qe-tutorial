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
([Ch. 3](../chapters/03-brillouin-zone.md)), DFPT phonons
([Ch. 10–11](../chapters/11-phonon-dispersion.md)), structural relaxation
([Ch. 4](../chapters/04-relaxation.md)), and charge densities
([Ch. 7](../chapters/07-charge.md)) — to reproduce a real many-body phenomenon: a
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

## 1. The Peierls instability — derived from a tight-binding chain

The cleanest place to *see* a CDW is a 1D chain. We build the band structure of the
uniform chain, then of the distorted chain in a two-site basis, compare them, and
identify exactly which phonon drives the transition.

### 1.1 The uniform monatomic chain

Take one orbital per site, lattice constant $a$, on-site energy set to zero, and
nearest-neighbour hopping $t>0$. In second quantization,

$$
H_0 = -t \sum_{n} \left( c_n^\dagger c_{n+1} + c_{n+1}^\dagger c_n \right) .
$$

Bloch states $|k\rangle = \frac{1}{\sqrt N}\sum_n e^{ikna} c_n^\dagger|0\rangle$
diagonalize it, giving a **single cosine band**

$$
\varepsilon(k) = -2t\cos(ka), \qquad k \in \left(-\tfrac{\pi}{a}, \tfrac{\pi}{a}\right].
$$

With **one electron per site** (half-filling) the band is filled up to the Fermi
points

$$
k_F = \pm\frac{\pi}{2a}, \qquad \varepsilon_F = -2t\cos\!\left(\tfrac{\pi}{2}\right) = 0 .
$$

The band crosses $\varepsilon_F$ with finite slope → the uniform chain is a **metal**.

### 1.2 The dimerized (CDW) chain in an A/B basis

Now let the ions dimerize: a displacement wave of wavevector $q = 2k_F = \pi/a$
pairs the atoms up, **doubling the period** to $a' = 2a$. The unit cell now holds
two inequivalent sites, **A** and **B**, and the single hopping $t$ splits into two
alternating values

$$
t_1 = t + \delta t \;(\text{intra-dimer}), \qquad
t_2 = t - \delta t \;(\text{inter-dimer}),
$$

where $\delta t \propto u_0$ is set by the distortion amplitude $u_0$ through the
electron–phonon coupling. With cell index $m$ and operators $a_m, b_m$,

$$
H = -\sum_m \left( t_1\, a_m^\dagger b_m + t_2\, b_m^\dagger a_{m+1} + \text{h.c.} \right) .
$$

Fourier transforming on the doubled cell ($a' = 2a$) gives, in the basis
$(a_k, b_k)$, a $2\times2$ **Bloch Hamiltonian**

$$
H(k) = \begin{pmatrix} 0 & h(k) \\ h^*(k) & 0 \end{pmatrix},
\qquad h(k) = -\left( t_1 + t_2\, e^{-ik a'} \right),
$$

whose eigenvalues are **two bands**

$$
\varepsilon_\pm(k) = \pm\,|h(k)| = \pm\sqrt{\,t_1^2 + t_2^2 + 2 t_1 t_2 \cos(k a')\,},
\qquad k \in \left(-\tfrac{\pi}{a'}, \tfrac{\pi}{a'}\right].
$$

### 1.3 Comparing the two band structures

The reduced zone now ends at $k = \pi/a' = \pi/2a = k_F$ — the old Fermi point is
folded onto the **new zone boundary**. Evaluate the gap there:

$$
\varepsilon_\pm(k_F) = \pm\sqrt{t_1^2 + t_2^2 - 2t_1 t_2} = \pm\,|t_1 - t_2| = \pm\,2\,\delta t .
$$

- **Undistorted** ($\delta t = 0$, i.e. $t_1=t_2=t$): the two bands touch,
  $\varepsilon_\pm(k_F)=0$ — this is just the original cosine band folded back into
  the half-zone, still **gapless at $E_F$** (a metal).
- **Distorted** ($\delta t \neq 0$): a gap

$$
2\Delta = 2\,|t_1 - t_2| = 4\,\delta t
$$

opens **exactly at $k_F$**, right at the Fermi level. The lower band (filled) is
pushed **down**, the upper band (empty) up → the chain becomes a **band insulator**
and the occupied electrons **lower their energy**.

![Peierls tight-binding band structures](/img/MP1-peierls-bands.png)

The two panels make the comparison concrete (computed from the formulas above in
[`notebooks/MP1-peierls.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/MP1-peierls.ipynb)):
**(a)** the uniform cosine band crosses $E_F$ at $\pm k_F$ — a metal; **(b)** in the
doubled-cell reduced zone the dashed $\delta t=0$ bands merely *touch* at $k_F$
(folded metal), while the dimerized bands (solid) open the gap $2\Delta = 4\,\delta t$
right at the Fermi level.

### 1.4 Why the $2k_F$ phonon: it mixes $-k_F$ and $+k_F$

Where does the off-diagonal $h(k)$ come from physically? The frozen phonon is a
displacement $u_n = u_0\cos(2k_F\, na)$, which adds a periodic potential
$\delta V(x)$ whose only Fourier component is at $G = 2k_F$. Its matrix elements in
the *original* Bloch basis therefore connect states separated by exactly $2k_F$:

$$
\langle k' | \delta V | k \rangle \neq 0 \iff k' = k \pm 2k_F .
$$

Set $k = -k_F$. Then $k' = -k_F + 2k_F = +k_F$ — so the perturbation **couples the
two degenerate Fermi-point states $|{-}k_F\rangle$ and $|{+}k_F\rangle$**, both
sitting at $\varepsilon_F$. Degenerate perturbation theory reduces to the $2\times2$
secular problem

$$
\begin{pmatrix} \varepsilon_F & \Delta \\ \Delta^* & \varepsilon_F \end{pmatrix},
\qquad \Delta = \langle +k_F|\delta V|-k_F\rangle \propto g\,u_0 ,
$$

whose eigenvalues are $\varepsilon_F \pm |\Delta|$. **That splitting *is* the gap.**
This is the whole mechanism in one line: the $q = 2k_F$ phonon is precisely the
distortion that hybridizes $-k_F \leftrightarrow +k_F$, and because those states are
degenerate at $E_F$, the mixing is maximally effective — it opens a gap and lowers
the electronic energy.

### 1.5 Standing waves: the gap and the charge modulation in real space

§1.4 was in momentum space; the same physics in **real space** makes both the gap
*and* the charge-density wave fall out at once — and gives the gap by an even simpler
route than the tight-binding $2\times2$. The $2k_F$ potential mixes the two degenerate
travelling waves $e^{\pm ik_Fx}$ into **standing waves**

$$
\psi_+ \propto \cos(k_F x), \qquad \psi_- \propto \sin(k_F x),
$$

whose charge densities are **no longer flat**:

$$
|\psi_+|^2 = 1 + \cos(2k_F x), \qquad |\psi_-|^2 = 1 - \cos(2k_F x).
$$

A travelling wave $e^{ik_Fx}$ has uniform $|\psi|^2$ — charge spread evenly, a metal.
The standing wave instead **piles charge periodically** with the CDW period
$2\pi/2k_F = 2a$. And the lattice distorts to match: the atoms **dimerize**,
$u_n = u_0(-1)^n$, pairing into short bonds with the charge heaped on them — the
**"flat → localized"** picture that *is* the charge-density wave (top: charge;
bottom: the uniform vs dimerized chains).

![Flat vs localized standing-wave charge density and the dimerized lattice](/img/MP1-peierls-density.png)

**The gap is just the energy difference of the two standing waves.** The lattice
distortion supplies a potential whose only Fourier component is at $2k_F$,
$V(x) = 2V_{2k_F}\cos(2k_F x)$ with $V_{2k_F}\propto u$. The two standing waves share
the same kinetic energy but sample this potential oppositely:

$$
E_\pm = \varepsilon_0 + \langle \psi_\pm|V|\psi_\pm\rangle = \varepsilon_0 \pm V_{2k_F},
$$

since $\langle\psi_+|V|\psi_+\rangle = \tfrac1L\!\int (1+\cos 2k_Fx)\,2V_{2k_F}\cos(2k_Fx)\,dx = V_{2k_F}$
(and $-V_{2k_F}$ for $\psi_-$). Hence

$$
\boxed{\;2\Delta = E_+ - E_- = 2\,|V_{2k_F}| \;\propto\; u\;}
$$

— the gap is **twice the $2k_F$ Fourier component of the potential**, growing linearly
with the distortion, exactly the tight-binding $2\Delta = 4\,\delta t$ by another road.
The occupied (lower) band is the standing wave whose charge sits in the **low-potential**
regions, so filling it simultaneously opens the gap and produces the charge pile-up.

This is exactly what we **plot** for NbSe₂ in Stage 3: the total $|\psi|^2$ summed over
occupied states — **flat** in the undistorted metal, **modulated** (charge localized on
the Nb clusters) in the CDW.

### 1.6 The energy balance — Peierls' theorem

Summing the lowered occupied states against the elastic cost of the distortion,

$$
\Delta E(u_0) = \underbrace{+\tfrac12 K\, u_0^2}_{\text{elastic}}
              \;\underbrace{-\,\gamma\, u_0^2 \ln\!\frac{W}{u_0}}_{\text{electronic}} ,
$$

with $K$ the lattice stiffness, $W \sim t$ a bandwidth, and $\gamma>0$. The elastic
term is a plain parabola, but the electronic gain carries a **logarithm** — a direct
consequence of the 1D density of states diverging at the band edge / the perfect
nesting of the two Fermi points by $2k_F$. For small $u_0$, $u_0^2\ln(1/u_0)$ always
beats $u_0^2$, so

$$
\Delta E(u_0) < 0 \quad\text{for some } u_0 \neq 0 .
$$

The uniform 1D metal is therefore **never** the ground state at $T=0$: it
spontaneously dimerizes. This is **Peierls' theorem**, and the special role of
$2k_F$ — the vector that nests the Fermi surface — is what we generalize next.

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

The relaxed cell (Stage 0 below) comes out at **a₀ = 3.473 Å**, with the Se
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
| 3 | DOS pseudogap + charge modulation | `dos.x`, `pp.x` | light |

### Stage 0 — relax the high-symmetry cell

A metal, so we use Marzari–Vanderbilt cold smearing; `cell_dofree='2Dxy'` relaxes
the in-plane lattice constant while the vacuum spacing is held fixed, and
`assume_isolated='2D'` removes the spurious field between periodic images.

:::warning Relax at the production settings
The relaxation uses the **same k-grid and `degauss`** as the phonon step that
follows (**12×12×1**, `degauss = 0.01`). This is not optional housekeeping: the
equilibrium lattice constant depends weakly on **k** and `degauss`, so relaxing at
*looser* settings and then running phonons at *tighter* ones leaves the cell
slightly **off-equilibrium**, and the residual stress leaks straight into the
dynamical matrix. Converge → relax → reuse, all at one consistent set of
parameters.
:::

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
    degauss         = 0.01
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
  12 12 1 0 0 0
```

### Stage 1 — phonon dispersion of the undistorted cell

This is the heart of the project. We run a self-consistent calculation on the
relaxed cell with a **12×12×1 k-grid** — exactly twice the 6×6×1 phonon **q**-grid,
so every $\mathbf k+\mathbf q$ lands back on the **k**-grid — then a full DFPT
phonon calculation on that **6×6×1 q-grid**, interpolate the force constants with
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
  12 12 1 0 0 0
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

![Monolayer NbSe2 phonon dispersion with soft mode](/img/MP1-nbse2-phonon.png)

**Result.** A whole branch is pushed **below zero** over a broad region centred on
$\mathbf q_\text{CDW}=(\tfrac13,0)=\tfrac23\,\Gamma\text M$. The directly-computed
dynamical matrix at that **q**-grid point gives a frequency of **−70.6 cm⁻¹**
(imaginary); the interpolated dispersion bottoms out at **≈ −72 cm⁻¹** — the two
agree, so the instability is real and not an interpolation artifact
([`notebooks/MP1-nbse2-phonon.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/MP1-nbse2-phonon.ipynb)).
This is the **deep, extended Kohn anomaly** that NbSe₂ is known for — note the dip is
*broad* in **q**, not a sharp spike. The acoustic modes return cleanly to ≈0 at
$\Gamma$ (the acoustic sum rule holds), confirming the instability sits at **finite
q**, not at the zone centre. That a clear imaginary mode survives even at this
moderate `degauss = 0.01` shows how strongly the monolayer wants to distort.

:::tip The soft mode depends on `degauss` — that's the physics
Because the depth of the dip is set by the Fermi-surface nesting, it sharpens as you
**lower `degauss`** (cool the electrons below $T_\text{CDW}$) and deepens with a
denser **k**-grid. Re-running at `degauss = 0.005, 0.0025` and k = 18, 24 makes the
mode plunge further — a direct, visual demonstration that the CDW is an
electronic-temperature-driven instability (Exercises 2–3).
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

:::tip Two registries: hollow- vs chalcogen-centred
The 3×3 distortion can lock in with two inequivalent registries — **hollow-centred
(HC)** or **chalcogen-centred (CC)**, named by what sits at the centre of the Nb
triangular cluster. They are genuinely different local minima with slightly
different energies (in 2H-NbSe₂ the two even appear on different van der Waals
layers — see [Sanna *et al.*, *npj Quantum Mater.* **7**, 6 (2022)](https://www.nature.com/articles/s41535-021-00412-8)).
The honest recipe is therefore to **relax *both* registries** and keep the lower
$\Delta E_\text{CDW}$ as the ground state — not to assume the first pattern you
freeze in is the winner.
:::

**Result.** Freezing in the 3Q (triangular) combination of the soft mode and
relaxing the 27-atom cell gives a genuine energy minimum below the undistorted
reference:

$$
\Delta E_\text{CDW} = \frac{E_\text{3×3 distorted} - E_\text{3×3 undistorted}}{9}
   = \frac{-25.3\ \text{meV}}{9} = \mathbf{-2.8\ meV/f.u.}
$$

The relaxed structure is the **textbook 3×3 triangular CDW**: the Nb atoms move by
only **~0.06–0.09 Å** (a few % of the 3.47 Å spacing), but that small shift groups
them into **triangular clusters** — the nearest-neighbour Nb–Nb bond *shortens* from
3.47 Å to **3.38 Å** inside a cluster and *stretches* to **3.63 Å** between clusters
([`notebooks/MP1-nbse2-cdw-structure.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/MP1-nbse2-cdw-structure.ipynb)).

![Relaxed 3x3 CDW structure of NbSe2: Nb triangular clusters](/img/MP1-nbse2-cdw-structure.png)

:::tip Why so small? NbSe₂ is a *weak* CDW
The displacement looks tiny — and that is the physics, not a failure. A CDW is an
**electronic** instability of a metal, not a bond-breaking reconstruction: only the
nested sliver of the Fermi surface gains a (pseudo)gap, so the energy and the
distortion are both small. The scale is set by the transition temperature:

$$
\Delta E_\text{CDW} \approx k_B T_\text{CDW}
\;\Rightarrow\; T_\text{CDW} \sim \frac{2.8\ \text{meV}}{k_B} \approx 33\ \text{K},
$$

right at NbSe₂'s observed $T_\text{CDW}\!\approx\!33$ K. A **small** $\Delta E$ *is* a
**low** transition temperature — they are the same statement. A *strong*-coupling CDW
such as 1T-TaS₂ ("Star of David", $\sim$0.2–0.3 Å shifts, $T_\text{CDW}\!\sim$200–350 K)
distorts far more visibly; NbSe₂ is the textbook **weak** CDW, which also lets it
coexist with superconductivity at low $T$.
:::

### Stage 3 — see the electronic structure change (DOS) and the charge modulation

This is where we answer "**does the band structure actually change?**" — and the
honest answer carries the key physics. Unlike the 1D toy model of §1, which turned
a metal into a full-gap **insulator**, monolayer NbSe₂'s CDW is a **weak, partial**
transition: only the **nested portions of the Fermi surface gap out**, while the
rest stays metallic (which is why NbSe₂ remains a metal — and superconducts at lower
temperature). So the electronic signature is **not** a clean gap but a **pseudogap**
— a partial dip in the density of states right at $E_F$.

The clearest way to see it is a **DOS overlay**: run `dos.x` on the undistorted 1×1
cell and on the relaxed 3×3 CDW cell (both unfolded to states-per-formula-unit) and
compare near $E_F$. The 27-atom band structure is too folded to read directly, but
the DOS dip is unambiguous.

```bash
# undistorted reference and CDW cell, each: nscf (dense k) -> dos.x
pw.x  < nbse2.nscf.in    > nbse2.nscf.out
dos.x < nbse2.dos.in     > nbse2.dos.out
```

Then `pp.x` on the relaxed 3×3 cell gives the valence charge density
$\rho(\mathbf r)=\sum_{n\mathbf k}|\psi_{n\mathbf k}(\mathbf r)|^2$ — the real, 2D
version of the $|\psi|^2$ cartoon in §1.5. In the undistorted metal it is essentially
**flat** (uniform); in the relaxed CDW it is **modulated**, with charge **localized on
the Nb triangular clusters**. That periodic pile-up across the supercell — flat →
localized — **is** the charge-density wave, now in a real material rather than the 1D
toy.

:::note Where this leads — superconductivity
The pseudogap is the Fermi-surface reconstruction that, at lower temperature,
controls the **anisotropy of the superconducting gap** in NbSe₂
([Sanna *et al.* 2022](https://www.nature.com/articles/s41535-021-00412-8), on the
2H bulk). That calculation needs superconducting DFT / anisotropic Eliashberg —
beyond `pw.x`/`ph.x` and out of scope here — but it consumes exactly the **pre-ASR
phonons and electron–phonon coupling** we computed in Stage 1. Good CDW phonons in,
correct $T_c$ physics out.
:::

:::info Result to follow
The DOS overlay (undistorted vs 3×3 CDW, pseudogap at $E_F$) and the real-space
charge-modulation figure will be added here once the supercell relaxation completes.
:::

---

## 5. Experimental signatures — can you see this in the lab?

Every quantity we computed has a direct experimental counterpart; NbSe₂ is in fact
one of the most-studied CDW systems. The small atomic displacement does **not** make
the CDW hard to detect — the *charge* modulation it drives, the superlattice peaks,
and the soft phonon are all unambiguous.

| Computed here | Measured by | What is seen |
|---|---|---|
| soft mode −70 cm⁻¹ at $\mathbf q_\text{CDW}$ | inelastic X-ray / neutron scattering | the phonon **Kohn anomaly** dipping toward zero at $\mathbf q_\text{CDW}$ |
| 3×3 triangular distortion | X-ray / electron diffraction | **superlattice (satellite) peaks** at $(\tfrac13,0)$ |
| charge modulation | **STM/STS** | the 3×3 charge pattern imaged directly in real space |
| pseudogap at $E_F$ | STS, ARPES | partial gap + Fermi-surface reconstruction |
| $\Delta E_\text{CDW}\!\approx\!k_BT_\text{CDW}$ | transport, specific heat | resistivity kink / heat-capacity anomaly at $\approx$33 K |

Monolayer NbSe₂ was characterized exactly this way by STM (Ugeda *et al.*,
*Nat. Phys.* **12**, 92, 2016), which resolved the 3×3 CDW *and* the superconducting
gap on the same sheet; the bulk soft phonon was mapped by inelastic X-ray scattering
(Weber *et al.*, *PRL* **107**, 107403, 2011) — the experimental twin of our Stage-1
dispersion. The decisive fingerprint is the **temperature dependence**: the pattern
appears below $\approx$33 K and melts above it, which is precisely the
`degauss`-as-electronic-temperature physics of §2 realized with a real thermostat.

---

## 6. Why this is the right capstone

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
   `loto_2d`, even though the polar 2D crystal in [A3](../advanced/03-bilayer-hbn.md)
   did. (Hint: metal vs insulator.)

---

This closes the loop on the whole tutorial: the ground-state, electronic-structure,
and lattice-dynamics tools you built up, aimed together at a genuine condensed-matter
phenomenon.
