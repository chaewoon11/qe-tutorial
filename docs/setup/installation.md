---
id: installation
title: Installation & Setup
sidebar_label: Installation
sidebar_position: 1
---

# Installation & Setup

What you need before Chapter 0.

## 1. Quantum ESPRESSO

You need QE ≥ 7.x with at least `pw.x`. Later chapters also use `dos.x`,
`projwfc.x`, `pp.x`, `bands.x`, `ph.x`, `q2r.x`, `matdyn.x`.

### Option A — module (HPC / managed workstation)

Many clusters ship QE as a module:

```bash
module avail 2>&1 | grep -i espresso     # find the module name
module load qe/7.5                        # or whatever it is called
which pw.x
```

### Option B — install yourself

Download from the [official site](https://www.quantum-espresso.org/) (or
`git clone` the [GitLab repo](https://gitlab.com/QEF/q-e)) and build:

```bash
./configure
make pw pp ph
```

See the [official install guide](https://www.quantum-espresso.org/Doc/user_guide/)
for MPI / parallel builds.

Verify the binary works:

```bash
echo "" | pw.x 2>&1 | head    # prints the PWSCF banner, then waits for input
```

## 2. Pseudopotentials

Every element needs a pseudopotential file (`.UPF`). This tutorial bundles the
GaAs pair so Chapter 0 runs out of the box:

```
code/pseudos/Ga.pbe-dn-kjpaw_psl.0.2.upf    # PBE PAW gallium (3d in valence)
code/pseudos/As.pbe-n-kjpaw_psl.0.2.upf     # PBE PAW arsenic
```

Both are PBE PAW potentials from QE's **PSlibrary** — a consistent pair (same
functional, same generation scheme), which matters: never mix pseudopotentials
made with different functionals.

For your own work, get a curated library:

- **[SSSP](https://www.materialscloud.org/discover/sssp)** — the standard
  accuracy/efficiency-tested set (recommended).
- **[PseudoDojo](http://www.pseudo-dojo.org/)** — norm-conserving sets with
  tested cutoffs.
- **PSlibrary** — ships with QE.

Drop the `.UPF` files into a directory and point `pseudo_dir` in your input
at it.

## 3. Visualization & plotting

- **[XCrySDen](http://www.xcrysden.org/)** or
  **[VESTA](https://jp-minerals.org/vesta/)** — structures, charge densities,
  Fermi surfaces (Chapter 7).
- **gnuplot** or **Python + matplotlib** — bands, DOS, convergence curves. The
  plotting notebooks in this tutorial use Python.

## 4. Get the tutorial files

Clone the repository so you have the runnable inputs and reference outputs
alongside you:

```bash
git clone https://github.com/chaewoon11/qe-tutorial.git
cd qe-tutorial/code
```

:::warning Where to run
`pw.x` calculations get expensive fast. The tiny 2-atom examples in early
chapters run in seconds on a laptop, but realistic cells, dense k-grids, and
especially DFPT phonons (Chapters 10–11) belong on an **HPC cluster**. Do not
run heavy jobs on a shared workstation.
:::
