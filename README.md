# Quantum ESPRESSO: Theory + Hands-On Tutorial

A from-scratch, in-depth tutorial for **Quantum ESPRESSO** (`pw.x` and friends).
Every chapter pairs the **theory** with a **hands-on calculation** you actually
run — theory → input-file dissection → run command → expected output → result
analysis → exercises.

> **Who this is for:** graduate students and researchers starting with
> plane-wave DFT who want to understand *why* each input flag is there, not just
> copy-paste. Some solid-state physics and command-line familiarity is assumed.

## What you'll learn

By the end you can independently run and interpret: ground-state SCF,
convergence studies, structural relaxation, band structures, DOS/PDOS, charge
densities, exchange-correlation choices (incl. DFT+U), magnetism, and phonons
via DFPT — and run them at scale on an HPC cluster.

## Curriculum

| # | Chapter | You'll learn | You'll run |
|---|---|---|---|
| 0 | What is QE / DFT in a nutshell | KS-DFT concept, the QE file ecosystem | a first Si SCF |
| 1 | Plane-wave DFT & the SCF loop | Bloch's theorem, plane-wave basis, self-consistency | `pw.x` scf, dissected |
| 2 | Pseudopotentials & cutoffs | NC vs US vs PAW, `ecutwfc`/`ecutrho`, SSSP | pseudopotential comparison |
| 3 | BZ sampling & convergence | Monkhorst–Pack grids, smearing | ecut/k convergence sweep |
| 4 | Structural optimization | Hellmann–Feynman forces, stress, BFGS | `relax` / `vc-relax` |
| 5 | Band structure | band folding, direct/indirect gaps, k-paths | `nscf` → `bands.x` |
| 6 | DOS & PDOS | density of states, orbital character, Löwdin charges | `dos.x`, `projwfc.x` |
| 7 | Charge density & visualization | charge density, ELF, STM | `pp.x` + XCrySDen/VESTA |
| 8 | Exchange-correlation functionals | LDA/GGA/metaGGA, DFT+U, hybrids, vdW | DFT+U on an oxide |
| 9 | Spin & magnetism | collinear & noncollinear magnetism, SOC | a magnetic system |
| 10 | DFPT & Γ phonons | linear response, the 2n+1 theorem, dynamical matrix | `ph.x` at Γ |
| 11 | Phonon dispersion | q-grids, IFCs, ASR, LO-TO splitting | `ph.x` → `q2r.x` → `matdyn.x` |
| 12 | Running QE at scale | parallelization levels, restarts, HPC job scripts | cluster submission |

Chapters 0–7 cover the ground state and electronic structure; 8–12 cover
functionals, magnetism, lattice dynamics, and HPC.

## How to use this repo

1. Start with [`setup/`](setup/) to install Quantum ESPRESSO and download
   pseudopotentials.
2. Work through `chapters/` in order. Each chapter is self-contained:
   ```
   chapters/NN-topic/
   ├── README.md     # the lesson (theory + walkthrough)
   ├── inputs/       # ready-to-run input files
   ├── ref/          # reference outputs and figures
   └── solutions/    # exercise solutions
   ```
3. Run the `inputs/`, compare against `ref/`, then try the exercises.

## Prerequisites

- Quantum ESPRESSO ≥ 7.x (`pw.x`, `ph.x`, `dos.x`, `projwfc.x`, `pp.x`, `bands.x`,
  `q2r.x`, `matdyn.x`) — see [`setup/`](setup/).
- A pseudopotential library (we use [SSSP](https://www.materialscloud.org/discover/sssp)).
- A visualizer ([XCrySDen](http://www.xcrysden.org/) or
  [VESTA](https://jp-minerals.org/vesta/)) and a plotting tool (gnuplot or
  Python/matplotlib).
- Heavy calculations should run on an HPC cluster — see [`docs/`](docs/) for
  Nurion / Stampede3 / Perlmutter guides.

## License & attribution

Tutorial material for educational use. Quantum ESPRESSO is distributed under the
GNU GPL; see the [official site](https://www.quantum-espresso.org/) and
[documentation](https://www.quantum-espresso.org/documentation/).
</content>
