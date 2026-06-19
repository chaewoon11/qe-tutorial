# Setup

What you need before Chapter 0.

## 1. Quantum ESPRESSO

You need QE ≥ 7.x with at least `pw.x`. Later chapters also use `dos.x`,
`projwfc.x`, `pp.x`, `bands.x`, `ph.x`, `q2r.x`, `matdyn.x`.

**Option A — module (HPC / managed workstation).** Many clusters ship QE:

```bash
module avail 2>&1 | grep -i espresso     # find the module name
module load qe/7.5                        # or whatever it is called
which pw.x
```

**Option B — install yourself.** Download from
<https://www.quantum-espresso.org/> (or `git clone` the
[GitLab repo](https://gitlab.com/QEF/q-e)) and build:

```bash
./configure
make pw pp ph
```

See the [official install guide](https://www.quantum-espresso.org/Doc/user_guide/)
for MPI/parallel builds.

Verify:

```bash
echo "test" | pw.x 2>&1 | head    # should print the PWSCF banner, then ask for input
```

## 2. Pseudopotentials

Every element needs a pseudopotential file (`.UPF`). This repo bundles **one**
so Chapter 0 runs out of the box:

```
setup/pseudos/Si_ONCV_PBE-1.1.upf     # norm-conserving PBE silicon (PseudoDojo)
```

For your own work, get a curated library:

- **[SSSP](https://www.materialscloud.org/discover/sssp)** — the standard
  accuracy/efficiency-tested set (recommended).
- **[PseudoDojo](http://www.pseudo-dojo.org/)** — norm-conserving sets with
  tested cutoffs.
- **PSlibrary** — ships with QE.

Drop the `.UPF` files into `setup/pseudos/` (or any directory) and point
`pseudo_dir` in your input at it.

## 3. Visualization & plotting

- **[XCrySDen](http://www.xcrysden.org/)** or **[VESTA](https://jp-minerals.org/vesta/)**
  — structures, charge densities, Fermi surfaces (Chapter 7).
- **gnuplot** or **Python + matplotlib** — bands, DOS, convergence curves.

## A note on where to run

`pw.x` calculations get expensive fast. The tiny 2-atom examples in early
chapters run in seconds on a laptop/workstation, but realistic cells, dense
k-grids, and especially DFPT phonons (Chapters 10–11) belong on an **HPC
cluster**. See [`../docs/`](../docs/) for Nurion / Stampede3 / Perlmutter job-script
guides. Do not run heavy jobs on a shared workstation.
