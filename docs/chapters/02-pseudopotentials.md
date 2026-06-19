---
title: "Chapter 2 — Pseudopotentials & cutoffs"
sidebar_label: "2 · Pseudopotentials & cutoffs"
sidebar_position: 3
---

# Chapter 2 — Pseudopotentials & cutoffs

Chapter 1 ended on a cliffhanger: plane waves can't resolve the wild
oscillations of wavefunctions near a nucleus without an astronomical cutoff, so
we replace the core with a **pseudopotential**. This chapter is about what that
replacement actually is, the three flavours you'll meet (NC / US / PAW), and how
the pseudopotential — together with a quick convergence study — fixes the two
cutoffs `ecutwfc` and `ecutrho`.

**Prerequisites:** [Chapter 1](./01-plane-wave-dft.md).

---

## 1. The frozen core, and why we "pseudize"

Two facts make core electrons a nuisance for a plane-wave code:

1. The valence wavefunctions oscillate rapidly near the nucleus (they must stay
   orthogonal to the tightly-bound core states). Representing those wiggles in
   plane waves would need cutoffs of thousands of Ry.
2. Core electrons barely change between, say, an isolated atom and the same atom
   in a crystal — chemistry is a valence affair.

So we make the **frozen-core approximation**: freeze the core electrons into the
nucleus and build an effective **pseudopotential** that (a) reproduces the true
valence wavefunctions *outside* a cutoff radius $r_c$, and (b) is smooth
*inside* it. The valence pseudo-wavefunctions are then nodeless and soft, and a
modest `ecutwfc` suffices. The price is that results are only as good as the
pseudopotential, and you must use one **consistent set** (same exchange–
correlation functional, ideally same generation scheme) across all elements.

:::warning Never mix functionals
A PBE calculation needs PBE pseudopotentials for *every* species. Mixing an LDA
pseudo with a PBE one (or two unrelated libraries) silently corrupts the result.
Our GaAs pair is PBE for both Ga and As, from the same PSlibrary generation.
:::

## 2. Three families: NC, US, PAW

| Family | Idea | `ecutwfc` | `ecutrho` (rule of thumb) | Notes |
|---|---|---|---|---|
| **Norm-conserving (NC)** | pseudo-wfc carries the same charge as the true one inside $r_c$ | highest | $4\times$ ecutwfc | simplest, "hardest", robust; e.g. ONCV |
| **Ultrasoft (US)** | relax norm-conservation → much softer wfc, add an augmentation charge | lowest | $8\text{–}12\times$ | cheap wfc, but needs a finer density grid |
| **PAW** | like US, but keeps the info to reconstruct the *all-electron* density | low | $8\text{–}12\times$ | US-like cost, more accurate; our GaAs set |

The trade-off is always the same: softer wavefunctions (lower `ecutwfc`) come at
the cost of a more complex, higher-frequency charge density (higher `ecutrho`).
Norm-conserving potentials need no augmentation, so `ecutrho = 4×ecutwfc` is
automatic; ultrasoft and PAW add an augmentation density that lives on the dense
grid, so they want a larger ratio. This is exactly why Chapter 0 had to set
`ecutrho` explicitly.

## 3. Read the pseudopotential header

A `.UPF` file is text, and its header tells you most of what you need. Look at
our gallium potential:

```bash
grep -iE "functional|pseudo_type|z_valence|core_correction|Suggested" \
     code/pseudos/Ga.pbe-dn-kjpaw_psl.0.2.upf
```

```text
    Functional:  SLA  PW   PBX  PBC
    Suggested minimum cutoff for wavefunctions:  33. Ry
    Suggested minimum cutoff for charge density: 134. Ry
             pseudo_type="PAW"
             core_correction="T"
             z_valence="13"
```

Decode it:

- **`Functional: SLA PW PBX PBC`** → this is **PBE**. Must match every other
  pseudo in the run.
- **`pseudo_type="PAW"`**, **`core_correction="T"`** → a PAW potential with a
  nonlinear core correction.
- **`z_valence="13"`** → 13 valence electrons: Ga 3d¹⁰4s²4p¹. The 3d shell is
  treated as *valence* (the `d` in the filename), which matters for accuracy in
  III–V semiconductors. As (`z_valence="5"`, 4s²4p³) keeps only the sp shell.
- **Suggested cutoffs 33 / 134 Ry** — the generator's recommended *minimums*.
  As suggests 20 / 103. You take the **largest** over all species and then
  *converge upward* from there.

The filename itself is a summary: `Ga` `.pbe` (functional) `-dn` (d-states in
valence + nonlinear core correction) `-kjpaw` (Kresse–Joubert PAW) `_psl.0.2`
(PSlibrary v0.2).

## 4. Hands-on: converge the cutoffs

"Suggested minimum" is a starting point, not a guarantee. The real test is to
raise each cutoff until the quantity you care about stops changing. Sweep
`ecutwfc` (holding `ecutrho` fixed and generous), recording the total energy:

```bash
cd code/00-first-scf
for ec in 20 30 40 50 60 70 80; do
  sed -e "s/ecutwfc       = 50.0/ecutwfc       = ${ec}.0/" \
      -e "s/conv_thr      = 1.0d-20/conv_thr      = 1.0d-10/" inputs/gaas.scf.in > /tmp/w.in
  echo "$ec $(grep '^!' <(pw.x < /tmp/w.in) | awk '{print $5}')"
done
```

(The convergence study itself doesn't need a `1e-20` SCF — `1e-10` is far below
the cutoff-induced differences and much faster.) The collected data are in
[`code/02-pseudopotentials/`](https://github.com/chaewoon11/qe-tutorial/tree/master/code/02-pseudopotentials);
[`notebooks/02-cutoff-convergence.ipynb`](https://github.com/chaewoon11/qe-tutorial/blob/master/notebooks/02-cutoff-convergence.ipynb)
turns them into convergence curves in **meV/atom** — the unit that matters:

![Cutoff convergence of GaAs](/img/02-cutoff-convergence.png)

Two different behaviours, and both are physics, not an accident:

- **`ecutwfc` (left)** — the total energy is **variational**: adding plane waves
  can only lower it, so the curve descends monotonically toward the true value.
  It is steep to ~40 Ry, then flattens. Our `ecutwfc = 50` sits a few meV/atom
  above the reference — fine for this tutorial; push to ~70 Ry for sub-meV
  production accuracy.
- **`ecutrho` (right)** — the energy is **not** variational in the density-grid
  cutoff (it wiggles up and down) and is only weakly sensitive: it stays within
  ~1 meV/atom across the whole range. `ecutrho` matters much more for **forces
  and stresses** than for the total energy, so picking it generously is cheap
  insurance rather than a tight convergence target.

:::tip Converge the observable you care about
Total energy is the easiest thing to converge, but you should ultimately
converge the *quantity you report* — a lattice constant, a band gap, a phonon
frequency. Those can need different cutoffs (and k-grids, Chapter 3) than the
energy. When in doubt, converge the property, not just $E$.
:::

## 5. Where to get pseudopotentials

Don't hand-pick potentials one at a time — use a curated, tested library:

- **[SSSP](https://www.materialscloud.org/discover/sssp)** (Standard Solid-State
  Pseudopotentials) — a per-element *best choice* curated for accuracy or for
  efficiency, each with a recommended `ecutwfc`/`ecutrho`. The usual default.
- **[PseudoDojo](http://www.pseudo-dojo.org/)** — systematically tested
  norm-conserving (and PAW) sets with hardness/accuracy metrics.
- **PSlibrary** — the US/PAW generation tables that ship with QE; our GaAs pair
  comes from here.

All of these give you a consistent set plus suggested cutoffs, which you then
confirm with a convergence study like §4.

## Exercises

1. **Pick the cutoff.** From the convergence curve, what `ecutwfc` keeps the
   total energy within **1 meV/atom** of the reference? Within 10 meV/atom?
2. **Trust but verify the header.** The Ga pseudo suggests `ecutwfc = 33`. Does
   running at exactly 33 Ry get you within 1 meV/atom (read it off the curve)?
   What does that tell you about "suggested minimum" values?
3. **Why does Ga have 13 valence electrons?** Write out its electron
   configuration and identify which shells this pseudopotential keeps in the
   valence. What would you expect to go wrong with a Ga pseudo that froze the 3d
   into the core?
4. **`ecutrho` and forces.** The total energy barely moved with `ecutrho`. Why
   might a *force* or *stress* calculation still demand a high `ecutrho`?
   (Hint: think about what the dense grid represents and how forces are
   derivatives.)

<details>
<summary>Show solutions</summary>

**1. Pick the cutoff.** Reading the left panel: ~70 Ry brings the energy under
1 meV/atom of the reference, while ~45–50 Ry already clears 10 meV/atom. The
tutorial's 50 Ry is the ~"few meV/atom" compromise — fast enough to iterate,
accurate enough to learn from.

**2. Trust but verify.** At 33 Ry the energy is still tens of meV/atom above the
converged value — the "suggested minimum" guarantees the pseudo is *usable*, not
that any *property* is converged. Always converge upward from the suggestion for
the accuracy you need.

**3. Ga 13 electrons.** Ga is [Ar]3d¹⁰4s²4p¹; this potential keeps **3d¹⁰ + 4s²
+ 4p¹ = 13** in the valence. Freezing the 3d into the core (a 3-electron Ga
pseudo) is cheaper but less transferable: the shallow 3d overlaps the valence
and responds to bonding, so a "small-core" 13-electron pseudo gives noticeably
better lattice constants and band structures for GaAs.

**4. `ecutrho` and forces.** The dense grid represents the **charge density** (and,
for PAW/US, the augmentation density), which has sharper features than the
wavefunctions. Forces and stresses are *derivatives* of the energy with respect
to atomic positions / cell, and derivatives amplify short-wavelength noise: an
under-resolved density grid produces "egg-box" force errors that the total
energy hides. Hence a generous `ecutrho` even when $E$ looks converged.

</details>

---

**Next:** [Chapter 3 — Brillouin-zone sampling & smearing](./03-brillouin-zone.md):
the other half of convergence — how dense a k-grid you need, and how metals
differ from semiconductors.
