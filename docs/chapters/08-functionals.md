---
title: "Chapter 8 — Exchange–correlation functionals & DFT+U"
sidebar_label: "8 · XC functionals & DFT+U"
sidebar_position: 9
---

# Chapter 8 — Exchange–correlation functionals & DFT+U

Every result so far used **PBE**. But the exchange–correlation functional
$E_\text{xc}[n]$ is the one piece of DFT we *approximate* (Chapter 1), and the
choice matters. This chapter surveys the ladder of functionals, shows how LDA and
PBE bracket the GaAs lattice constant, and introduces **DFT+U** for the strongly
correlated transition-metal oxide NiO — where plain PBE qualitatively fails.

**Prerequisites:** [Chapter 4](./04-relaxation.md) (equation of state) and
[Chapter 6](./06-dos.md) (DOS); a glance ahead at Chapter 9 helps for the
magnetism.

---

## 1. Jacob's ladder of functionals

Approximations to $E_\text{xc}[n]$ are usually ranked by what local information
they use — "Jacob's ladder", climbing toward chemical accuracy at rising cost:

| Rung | Uses | Examples | Typical behaviour |
|---|---|---|---|
| **LDA** | density $n$ | PZ, PW | over-binds: lattice constants **too small**, gaps too small |
| **GGA** | $n$ and $\nabla n$ | **PBE**, PBEsol | PBE over-corrects: lattices **too large**; PBEsol aims for the middle |
| **meta-GGA** | + kinetic-energy density $\tau$ | SCAN, r²SCAN | better geometries & energetics, pricier |
| **hybrid** | + exact (Fock) exchange | HSE, PBE0 | much better **gaps**, expensive (Fock in plane waves) |

Orthogonal to the ladder are **dispersion corrections** for van der Waals
bonding (`vdw_corr = 'dft-d3'`, or the vdW-DF functionals), essential for layered
and molecular systems where GGA has no long-range attraction.

:::note No rung is universally "best"
Higher rungs are usually more accurate but not always, and they cost more. The
right functional depends on the property and the material — which is exactly why
you should know how sensitive your answer is to the choice.
:::

## 2. Hands-on: LDA vs PBE on GaAs

The cleanest demonstration is the equation of state. We repeat the Chapter 4 EOS
with **LDA** pseudopotentials (`Ga.pz-dn` + `As.pz-n`, the PZ/LDA members of the
same PSlibrary family as our PBE pair) and compare — same cell, cutoffs, and
k-grid, only the functional changes:

```bash
cd code/08-functionals
# LDA EOS: loop the LDA scf input over celldm, grep the total energy
for a in 10.0 10.2 10.4 10.5 10.6 10.7 10.8 11.0; do
  sed "s/celldm(1) *= .*/celldm(1) = $a/" gaas.lda.scf.in > /tmp/e.in
  echo "$a $(grep '^!' <(pw.x < /tmp/e.in) | tail -1 | awk '{print $5}')"
done
```

:::warning Match the pseudopotential to the functional
A proper functional comparison uses pseudopotentials *generated* for each
functional — LDA pseudos for an LDA run, PBE pseudos for a PBE run. Overriding
`input_dft` on a mismatched pseudo is a shortcut that can mislead. Here we use the
PZ and PBE versions of the same PSlibrary potentials.
:::

[`notebooks/08-lda-vs-pbe.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/08-lda-vs-pbe.ipynb)
fits both curves:

![GaAs equation of state: LDA vs PBE](/img/08-lda-pbe-eos.png)

The two minima straddle experiment:

| Functional | $a_0$ (Å) | error vs exp (5.653 Å) |
|---|---|---|
| **LDA** | 5.612 | −0.7 % (too small) |
| **PBE** | 5.756 | +1.8 % (too large) |

LDA over-binds and contracts the cell; PBE over-corrects and expands it. The
experimental value sits in between — a textbook illustration that the functional,
not the convergence, sets the structure.

## 3. Why PBE fails for NiO: strong correlation

GGA is fine for the sp-bonded GaAs. It breaks down for **transition-metal
oxides** with partially-filled, localized **d** shells. The culprit is the
**self-interaction error**: semilocal functionals let an electron spuriously
repel itself, which *over-delocalizes* the d electrons and washes out the strong
on-site Coulomb repulsion. The result is qualitatively wrong — PBE predicts NiO
to be a small-gap (almost metallic) weak antiferromagnet, when it is really a
wide-gap (~4 eV) Mott/charge-transfer insulator.

**DFT+U** patches this cheaply. It adds a Hubbard-like penalty that pushes the
d-shell occupations toward integer values (empty or filled), restoring the
localization:

$$
E_\text{DFT+U} = E_\text{DFT} + \frac{U}{2}\sum_{I,\sigma}
\mathrm{Tr}\!\left[\mathbf{n}^{I\sigma}\left(1-\mathbf{n}^{I\sigma}\right)\right],
$$

where $\mathbf{n}^{I\sigma}$ is the occupation matrix of the localized d orbitals
on atom $I$ and $U$ is the effective Hubbard parameter (it can be computed from
linear response / the HP code, or taken from the literature — here $U=5$ eV for
Ni 3d). In modern QE this is set with the **`HUBBARD` card**.

## 4. Hands-on: NiO with DFT+U

NiO is rocksalt, ordered **antiferromagnetically (type II)**: ferromagnetic
(111) planes of Ni alternate in spin. The magnetic cell holds two Ni (opposite
spins) and two O, so we use **two Ni species** (`Ni1`, `Ni2`) with opposite
`starting_magnetization` and `nspin = 2`:

```fortran title="code/08-functionals/nio.pbeu.scf.in"
&control
    calculation = 'scf'
    prefix      = 'nio'
    outdir      = './out'
    pseudo_dir  = '../pseudos'
    verbosity   = 'high'
/
&system
    ibrav = 0
    celldm(1) = 7.881
    nat = 4
    ntyp = 3
    ecutwfc = 75.0
    ecutrho = 600.0
    nbnd = 30
    nspin = 2
    starting_magnetization(1) =  0.6      ! Ni1 up
    starting_magnetization(2) = -0.6      ! Ni2 down
    starting_magnetization(3) =  0.0      ! O
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.01
/
&electrons
    mixing_beta = 0.3
    conv_thr = 1.0d-8
/
ATOMIC_SPECIES
  Ni1 58.6934 Ni.pbe-spn-kjpaw_psl.1.0.0.UPF
  Ni2 58.6934 Ni.pbe-spn-kjpaw_psl.1.0.0.UPF
  O   15.9994 O.pbe-n-kjpaw_psl.1.0.0.UPF

CELL_PARAMETERS alat
  0.50 0.50 1.00
  0.50 1.00 0.50
  1.00 0.50 0.50

ATOMIC_POSITIONS crystal
  Ni1 0.00 0.00 0.00
  Ni2 0.50 0.50 0.50
  O   0.25 0.25 0.25
  O   0.75 0.75 0.75

K_POINTS automatic
  4 4 4 0 0 0

HUBBARD {atomic}
U Ni1-3d 5.0
U Ni2-3d 5.0
```

Run it twice — once with `U` set to ~0 (plain PBE) and once at `U = 5.0` eV — and
compare the magnetic moment and the gap:

```bash
cd code/08-functionals
pw.x < nio.pbe.scf.in  > nio.pbe.out     # HUBBARD U = 0.0001 (≈ PBE)
pw.x < nio.pbeu.scf.in > nio.pbeu.out    # HUBBARD U = 5.0 eV
grep -E "absolute magnetization|highest occupied" nio.pbe.out nio.pbeu.out
```

The contrast is the whole point:

| | PBE (U≈0) | PBE+U (U=5 eV) | experiment |
|---|---|---|---|
| Ni magnetic moment | ~1.5 μB | ~1.7 μB | ~1.6–1.9 μB |
| band gap | 1.0 eV | 2.65 eV | ~4.0 eV |

(The moment per Ni is half the `absolute magnetization` QE reports — 3.08 →
3.48 μB/cell for two Ni.) Turning on `U` does exactly what the theory promised:
it **opens the gap** (1.0 → 2.65 eV) and **pushes the Ni moment toward
experiment** (1.5 → 1.7 μB), localizing the d electrons. At `U = 5 eV` the gap
still undershoots the ~4 eV experiment; a larger `U` (~7–8 eV, the
linear-response value for this system) gets closer — which is the point of
Exercise 2.

:::tip `U` is not a free fudge factor
A different `U` gives a different gap, so don't tune it to taste. Either take a
well-justified literature value for your material/pseudo, or compute it from
first principles (linear response, or QE's **HP** code). Report the `U` you used.
:::

## Exercises

1. **PBEsol prediction.** Given LDA gives 5.61 Å and PBE 5.76 Å, where would you
   expect **PBEsol** (a GGA retuned for solids) to land? Why do solid-state
   workflows often prefer it over PBE for geometries?
2. **The role of `U`.** Re-run NiO at `U = 0`, `3`, `5`, `8` eV. Plot the gap vs
   `U`. Is the dependence linear? What does that say about quoting a gap "from
   DFT+U"?
3. **AFM vs FM.** Set both Ni `starting_magnetization` to `+0.6` (ferromagnetic)
   and compare the total energy with the antiferromagnetic ground state. Which is
   lower, and by how much?
4. **Why two Ni species?** Explain why the input declares `Ni1` and `Ni2` as
   separate `ATOMIC_SPECIES` even though they are the same element. What breaks if
   you use a single `Ni` type?

<details>
<summary>Show solutions</summary>

**1. PBEsol.** It should land **between** LDA and PBE, near experiment (~5.65 Å):
PBEsol restores the LDA-like gradient behaviour for slowly-varying densities,
fixing PBE's over-expansion of solids. That is why it is a popular default for
lattice constants and phonons.

**2. The role of `U`.** The gap grows roughly **linearly** with `U` over a wide
range — which is exactly why `U` cannot be chosen by eye: any gap you want is a
`U` away. Quote the `U` and its justification (literature or computed), never a
gap "from DFT+U" alone.

**3. AFM vs FM.** The **antiferromagnetic** type-II ordering is the ground state;
the ferromagnetic arrangement is higher in energy (the energy difference maps to
the magnetic exchange couplings). Finding AFM lower is a sanity check that the
magnetic setup is physical.

**4. Two Ni species.** They must be distinct types so they can carry **opposite**
`starting_magnetization` and (optionally) independent occupation matrices /
`U`. With a single `Ni` type, every Ni would start with the same spin — you would
converge to the ferromagnetic (or non-magnetic) state and never find the AFM
ground state.

</details>

---

**Next:** [Chapter 9 — Spin & magnetism](./09-magnetism.md): collinear `nspin=2`
in full, noncollinear magnetism, and spin–orbit coupling.
