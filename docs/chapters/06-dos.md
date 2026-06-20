---
title: "Chapter 6 — Density of states & PDOS"
sidebar_label: "6 · DOS & PDOS"
sidebar_position: 7
---

# Chapter 6 — Density of states & PDOS

A band structure shows eigenvalues along a *path*; the **density of states**
(DOS) counts them over the *whole* Brillouin zone — how many states sit at each
energy. Projecting that count onto atomic orbitals (the **PDOS**) then tells you
*what each band is made of*. This chapter computes both for GaAs with `dos.x`
and `projwfc.x`, and connects the peaks back to the bands of Chapter 5.

**Prerequisites:** [Chapter 5](./05-bands.md).

---

## 1. What the DOS is

The density of states is the number of Kohn–Sham states per unit energy,

$$
g(E) = \sum_n \int_\text{BZ} \frac{d\mathbf{k}}{V_\text{BZ}}\,
\delta\!\left(E - \varepsilon_{n\mathbf{k}}\right).
$$

Integrating $g(E)$ up to the Fermi level gives the total number of electrons;
peaks in $g(E)$ are energies with many available states (flat bands), and a gap
shows up as a region where $g(E)=0$. Because it is a BZ integral of a sharp
$\delta$-function, the DOS needs a **dense k-grid** — much denser than the SCF —
and a good integration method (**tetrahedra**, which needs no smearing width).

## 2. Total DOS: scf → dense nscf → `dos.x`

```bash
cd code/06-dos
pw.x  < gaas.scf.in   > scf.out      # self-consistent density (8×8×8)
pw.x  < gaas.nscf.in  > nscf.out     # dense grid for integration (16×16×16)
dos.x < gaas.dos.in   > dos.out      # -> gaas.dos
```

The nscf input differs from the SCF in the ways that matter for a DOS:
`occupations = 'tetrahedra'`, more bands, and a denser uniform grid.

```fortran title="code/06-dos/gaas.nscf.in"
&control
    calculation = 'nscf'
    prefix      = 'gaas'
    outdir      = './out'
    pseudo_dir  = '../pseudos'
    verbosity   = 'high'
/
&system
    ibrav       = 0
    celldm(1)   = 10.6829
    nat         = 2
    ntyp        = 2
    ecutwfc     = 50.0
    ecutrho     = 400.0
    occupations = 'tetrahedra'   ! clean BZ integration, no smearing parameter
    nbnd        = 24             ! cover well into the conduction bands
/
&electrons
    conv_thr    = 1.0d-10
    mixing_beta = 0.7
/
ATOMIC_SPECIES
  Ga  69.723   Ga.pbe-dn-kjpaw_psl.0.2.upf
  As  74.9216  As.pbe-n-kjpaw_psl.0.2.upf

CELL_PARAMETERS alat
  -0.50   0.00   0.50
   0.00   0.50   0.50
  -0.50   0.50   0.00

ATOMIC_POSITIONS alat
  Ga  0.00  0.00  0.00
  As  0.25  0.25  0.25

K_POINTS automatic
  16 16 16 0 0 0                 ! dense, uniform (not a path!)
```

`dos.x` then bins the eigenvalues into `gaas.dos` (columns: energy, $g(E)$,
integrated DOS), and reports the Fermi level in the header.

## 3. Projected DOS: `projwfc.x`

`projwfc.x` projects each Kohn–Sham state onto the **atomic orbitals** of the
pseudopotentials (a Löwdin projection) and resolves the DOS by atom and orbital:

```bash
projwfc.x < gaas.proj.in > proj.out
```

It writes one file per atom-and-orbital, e.g.
`gaas.pdos_atm#1(Ga)_wfc#3(d)`, plus a `gaas.pdos_tot`. It also prints **Löwdin
charges** — each state's projection integrated over the occupied manifold, i.e.
how many electrons of each orbital character sit on each atom.

## 4. Hands-on: the GaAs DOS, decomposed

[`notebooks/06-dos.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/06-dos.ipynb)
plots the total DOS with the per-element, per-orbital projections, aligned to the
valence-band maximum:

![GaAs density of states (PBE)](/img/06-dos.png)

Every feature maps onto a band from Chapter 5:

- The narrow spike near **−14.8 eV** is **Ga 3d** — sharp because the semicore
  3d states are flat bands (almost no dispersion → almost all states at one
  energy).
- The structure around **−12 to −7 eV** is the **As 4s** (and some Ga s) lower
  valence.
- The broad block from **~−7 eV up to the VBM** is dominated by **As 4p** — the
  bonding states that form the top of the valence band.
- Above the **gap** (the zero-DOS strip just above 0) the conduction band is a
  mix of Ga and As **s/p** character.

This is exactly the payoff of the PDOS: the *same* total DOS, but now you can see
that the valence-band top is As-p-like and the deep flat peak is Ga-d-like —
chemistry you cannot read off the bare band structure.

### Löwdin charges

`projwfc.x` also prints the integrated projections:

```text
Lowdin Charges:
  Atom #  1 (Ga): total = 12.96,  s = 1.08,  p = 1.88,  d = 9.99
  Atom #  2 (As): total =  5.01,  s = 1.44,  p = 3.57
Spilling Parameter:  0.0019
```

These confirm the orbital occupancies: the **Ga 3d shell is essentially full**
(d ≈ 9.99, as a semicore shell should be), while the bonding charge sits in the
As 4p (3.57) and the s/p manifolds. The tiny **spilling parameter (0.0019)** says
the atomic-orbital basis captures the occupied states almost perfectly — when
spilling is large, the projection (and the PDOS) is less trustworthy. (Löwdin
charges are projections, not a unique "ionic charge" — read them as trends, not
absolute oxidation states.)

## Exercises

1. **Where do the 18 electrons go?** Use the integrated DOS (third column of
   `gaas.dos`) at the VBM. Does it equal 18? Why must it?
2. **Tetrahedra vs smearing.** Re-run the nscf with `occupations='smearing'`,
   `degauss=0.01` and compare the DOS near the band edges. Which gives sharper
   edges, and why is tetrahedra preferred for an insulator DOS?
3. **Read the chemistry.** From the PDOS, which orbital dominates the *top* of
   the valence band, and which dominates the *bottom* of the conduction band?
   What does that imply about the states a hole vs. an electron would occupy?
4. **k-grid for DOS.** Halve the nscf grid to 8×8×8 and re-plot. Where does the
   DOS get noticeably worse, and why do sharp features (like the Ga 3d peak)
   suffer most?

<details>
<summary>Show solutions</summary>

**1. Where the electrons go.** The integrated DOS at the VBM equals **18** — the
number of valence electrons per cell (Ga 13 + As 5, doubly counted as 9 filled
bands × 2). It must, because integrating $g(E)$ to the highest occupied level
counts every occupied state.

**2. Tetrahedra vs smearing.** The tetrahedron method gives sharper, more faithful
band edges and van Hove singularities; Gaussian smearing artificially broadens
them by `degauss` and can smear states into the gap. For an insulator DOS,
tetrahedra (no width parameter) is the cleaner choice.

**3. Read the chemistry.** The valence-band top is **As 4p**; the conduction-band
bottom is **s-like** (mixed Ga/As s). So a *hole* lives in As p-character states,
while an injected *electron* occupies the s-like Γ conduction minimum — relevant
to doping and optical transitions.

**4. k-grid for DOS.** A coarser grid makes the DOS jagged (sampling noise),
worst where the DOS varies fastest — sharp peaks and band edges. The Ga 3d spike,
being nearly a delta function, is the most sensitive: too few k-points render it
as scattered spikes instead of one clean peak.

</details>

---

**Next:** [Chapter 7 — Charge density & visualization](./07-charge.md): `pp.x`,
real-space densities, and looking at bonds with XCrySDen / VESTA.
