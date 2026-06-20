---
id: intro
title: Quantum ESPRESSO — Theory + Hands-On
slug: /
sidebar_position: 1
sidebar_label: Introduction
---

# Quantum ESPRESSO: Theory + Hands-On

A from-scratch, in-depth tutorial for **Quantum ESPRESSO** (`pw.x` and friends).
Every chapter pairs the **theory** with a **hands-on calculation** you actually
run — theory → input-file dissection → run command → expected output → result
analysis → exercises.

:::note Who this is for
Graduate students and researchers starting with plane-wave DFT who want to
understand *why* each input flag is there, not just copy-paste. Some
solid-state physics and command-line familiarity is assumed.
:::

## What you'll learn

By the end you can independently run and interpret: ground-state SCF,
convergence studies, structural relaxation, band structures, DOS/PDOS, charge
densities, exchange–correlation choices (incl. DFT+U), magnetism, and phonons
via DFPT — and run them at scale on an HPC cluster.

## Curriculum

| # | Chapter | You'll learn | You'll run |
|---|---|---|---|
| 0 | What is QE / DFT in a nutshell | KS-DFT concept, the QE file ecosystem | a first Si SCF |
| 1 | Plane-wave DFT & the SCF loop | Bloch's theorem, plane-wave basis, self-consistency | `pw.x` scf, dissected |
| 2 | Pseudopotentials & cutoffs | NC vs US vs PAW, `ecutwfc`/`ecutrho` | pseudopotential comparison |
| 3 | BZ sampling & convergence | Monkhorst–Pack grids, smearing | ecut/k convergence sweep |
| 4 | Structural optimization | Hellmann–Feynman forces, stress, BFGS | `relax` / `vc-relax` |
| 5 | Band structure | band folding, direct/indirect gaps, k-paths | `nscf` → `bands.x` |
| 6 | DOS & PDOS | density of states, orbital character, Löwdin charges | `dos.x`, `projwfc.x` |
| 7 | Charge density & visualization | charge density, ELF, STM | `pp.x` + XCrySDen/VESTA |
| 8 | Exchange–correlation functionals | LDA/GGA/metaGGA, DFT+U, hybrids, vdW | DFT+U on an oxide |
| 9 | Spin & magnetism | collinear & noncollinear magnetism, SOC | a magnetic system |
| 10 | DFPT & Γ phonons | linear response, the 2n+1 theorem, dynamical matrix | `ph.x` at Γ |
| 11 | Phonon dispersion | q-grids, IFCs, ASR, LO-TO splitting | `ph.x` → `q2r.x` → `matdyn.x` |
| 12 | Running QE at scale | parallelization levels, restarts, HPC job scripts | cluster submission |

Chapters 0–7 cover the ground state and electronic structure; 8–12 cover
functionals, magnetism, lattice dynamics, and HPC.

:::info Status
**All 13 chapters (0–12) are complete and runnable** — the full curriculum, from
your first SCF to phonon dispersions and running QE at scale on HPC, every step
verified against real Quantum ESPRESSO output.
:::

## How to use this site

1. Work through [**Setup**](setup/installation.md) to install Quantum ESPRESSO and
   get pseudopotentials.
2. Read the chapters in order. Each is self-contained: theory, an annotated
   input file, the run command, the output explained line by line, and
   exercises.
3. The runnable input files, reference outputs, and the bundled GaAs
   pseudopotentials live in the
   [`code/`](https://github.com/chaewoon11/qe-tutorial/tree/master/code)
   directory of the repository — clone it and follow along.

## Prerequisites

- Quantum ESPRESSO ≥ 7.x (`pw.x`, `ph.x`, `dos.x`, `projwfc.x`, `pp.x`,
  `bands.x`, `q2r.x`, `matdyn.x`).
- A pseudopotential library — we use
  [SSSP](https://www.materialscloud.org/discover/sssp).
- A visualizer ([XCrySDen](http://www.xcrysden.org/) or
  [VESTA](https://jp-minerals.org/vesta/)) and a plotting tool
  (gnuplot or Python + matplotlib).

Quantum ESPRESSO is distributed under the GNU GPL; see the
[official site](https://www.quantum-espresso.org/).
