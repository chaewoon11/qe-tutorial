---
title: "Chapter 12 — Running QE at scale"
sidebar_label: "12 · Running QE at scale"
sidebar_position: 13
---

# Chapter 12 — Running QE at scale

The phonon jobs of Chapters 10–11 were the first that genuinely needed a cluster.
This final chapter is the practical layer: how `pw.x` (and friends) parallelize
over **several independent levels**, how to choose them, how to **restart** a job
that hits the wall-clock limit, and what a real HPC **job script** looks like.

**Prerequisites:** the rest of the tutorial — this is about running the
calculations you now understand, faster.

---

## 1. QE's parallelization levels

QE doesn't have one "number of cores" knob; it splits the work along a hierarchy
of (mostly) independent levels. From coarsest (least communication) to finest
(most):

| Level | Flag | Splits | Communication | Best for |
|---|---|---|---|---|
| **images** | `-nimage` | independent images (NEB, ph.x q/irreps) | almost none | NEB, phonon work |
| **pools** | `-npool` | the k-points | low | many k-points (metals, small cells, dense grids) |
| **bands** | `-nband` | groups of Kohn–Sham bands | medium | many bands (big systems) |
| **plane waves** | *(the rest)* | G-vectors / the FFT grid | **high** (FFTs) | the base level; reduces memory |
| **task groups** | `-ntg` | FFT task groups | high | many ranks per pool |
| **linear algebra** | `-ndiag` | the dense subspace diagonalization (ScaLAPACK) | high | large band counts |

On top of the MPI levels, QE is also **hybrid**: each MPI rank can run several
**OpenMP threads** (`OMP_NUM_THREADS`).

The total MPI ranks factorize across the levels. Within a pool, the ranks left
after pooling do the plane-wave/FFT work (and, optionally, the band and
linear-algebra splits). The art is putting cores where the communication is
cheapest.

## 2. Choosing the flags

A few rules carry most of the way:

- **`-npool` is your first lever** when you have many k-points. K-points are
  independent, so pools scale almost linearly with little communication.
  Constraints: `npool` must **divide the total number of MPI ranks**, and should
  be **≤ the number of (irreducible) k-points**.
- The ranks *within* each pool (`ranks_per_pool = total_ranks / npool`) split the
  **plane waves / FFT** — the communication-heavy part. Don't make pools so small
  that the FFT can't be distributed, but don't over-parallelize PWs across nodes
  either (FFT all-to-all is expensive over the network).
- **`-ndiag`** (the ScaLAPACK grid for subspace diagonalization) must be a
  **perfect square ≤ ranks_per_pool**. For small/medium systems `-ndiag 1`
  (serial diag) is often fastest; it matters for large band counts.
- **Memory:** pools **duplicate** data (each pool holds the whole G-vector set),
  so *more pools = more memory*. Plane-wave parallelization, in contrast,
  *distributes* the G-vectors and **reduces** per-rank memory. If you hit
  out-of-memory, use **fewer pools / more PW ranks** (or more nodes).

**Worked example (our Nurion phonon job).** One KNL node = 68 MPI ranks, run as
`mpirun -np 68 ph.x -npool 4`: 4 k-point pools × 17 ranks/pool on the FFT. That
matched the job's k-points and kept each pool large enough for the plane-wave
work.

## 3. Hands-on: how far does a small job scale?

A scaling test on the light GaAs SCF (`code/12-hpc/gaas.scf.in`, 16 irreducible
k-points), spreading the k-points across cores with `-npool = cores`:

```bash
export OMP_NUM_THREADS=1
for n in 1 2 4 8 16; do
  mpirun -np $n pw.x -npool $n < gaas.scf.in > scf_$n.out
done
```

| cores (`-npool`) | wall time | speedup | efficiency |
|---|---|---|---|
| 1 | 8.0 s | 1.0× | 100 % |
| 2 | 5.7 s | 1.4× | 71 % |
| 4 | 3.7 s | 2.2× | 55 % |
| 8 | 2.9 s | 2.8× | 35 % |
| 16 | 2.8 s | 2.9× | 18 % |

![GaAs SCF strong scaling](/img/12-scaling.png)

The speedup is real but **saturates fast**: by ~8 cores this 2-atom cell is so
cheap that communication and serial overhead dominate, and 16 cores is barely
faster than 8 (efficiency down to 18 %). That is the flip side of "use HPC for
the big jobs" — **don't over-parallelize a small one**, or you waste most of the
cores. Match the core count to the work.

:::tip Scaling is system-specific — always test first
This cell is tiny and k-point-rich, so `-npool` carries it (16 independent
k-points spread across cores with little communication) but saturates early. A
*large* cell with *few* k-points scales the other way — through plane-wave/band
parallelism and OpenMP — and keeps scaling to many more cores. Never assume:
run a short scaling test on *your* system before requesting a big allocation.
:::

## 4. Restarts and the wall clock

HPC jobs have a wall-time limit; long runs must **checkpoint and resume**.

- **`pw.x`:** set `max_seconds` in `&control` a bit below your wall time so the
  SCF stops cleanly and writes the restart files; resume with
  `restart_mode = 'restart'`. The `outdir` (with `prefix.save`) must persist
  between submissions.
- **`ph.x`:** the DFPT run is resumable with **`recover = .true.`** — it picks up
  from the last finished representation/q-point. This is essential for big
  dispersions that exceed one wall-clock window.

:::warning Never wipe the scratch between resubmits
The single most common way to lose days of DFPT: a job script that starts with
`rm -rf *.save _ph0 save`. That deletes exactly the recover state, so a
walltime-killed `ph.x` restarts from scratch. Make the script **recover-safe**
(don't clear those directories); force a fresh run by clearing the directory *by
hand* when you actually mean to.
:::

## 5. A real job script

Every job is the same shape: scheduler directives, environment, then `mpirun`.
This is the **PBS** script we used for the Chapter 11 dispersion on Nurion (KNL):

```bash title="job.sh (PBS / Nurion)"
#!/bin/bash
#PBS -N gaas_phdisp
#PBS -q debug
#PBS -l select=1:ncpus=68:mpiprocs=68:ompthreads=1
#PBS -A <allocation>
#PBS -l walltime=06:00:00
#PBS -j oe
#PBS -o run.out

cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi_21.2 impi/oneapi_21.2 craype-mic-knl
export QE=$HOME/qe-7.5-knl/q-e-qe-7.5/bin

mpirun -np 68 $QE/pw.x -npool 4 < gaas.scf.in     > gaas.scf.out
mpirun -np 68 $QE/ph.x -npool 4 < gaas.ph.disp.in > gaas.ph.disp.out
```

A **SLURM** cluster (e.g. Perlmutter / Stampede3) is the same idea with
`#SBATCH` directives and `srun` instead of `mpirun`. Our lab's cluster-specific
guides (modules, queues, gotchas) are in
[`reference/hpc-guides/`](https://github.com/chaewoon11/qe-tutorial/tree/master/reference/hpc-guides)
(Nurion / Stampede3 / Perlmutter).

:::note Where to run what
A rough guide from this tutorial: SCF / relax / bands / DOS on a workstation or a
few cores; **DFPT phonons, dense grids, and large cells on HPC**. The 2-atom Γ
phonon (Chapter 10) was minutes on one node; a full dispersion (Chapter 11) is
many such jobs — exactly the kind of thing the scheduler is for.
:::

## Exercises

1. **Pick a layout.** You have 96 MPI ranks and a calculation with 12
   irreducible k-points. Give a sensible `-npool`, the resulting ranks-per-pool,
   and a valid `-ndiag`. What changes if there were only 2 k-points?
2. **Out of memory.** A job dies with an allocation error. You're running
   `-npool 8`. Does *increasing* or *decreasing* `npool` help the memory, and
   why?
3. **Walltime hit.** Your `ph.x` dispersion is killed at the 12 h limit
   two-thirds done. What one input flag lets you finish it in a second job, and
   what must you *not* do in the resubmission script?
4. **The right machine.** For each — a 2-atom SCF, a 200-atom relaxation, a phonon
   dispersion — say whether you'd run on a workstation or HPC, and which
   parallel level you'd lean on.

<details>
<summary>Show solutions</summary>

**1. Pick a layout.** `-npool 12` (divides 96, equals the k-points) → 8
ranks/pool for the FFT; `-ndiag` a perfect square ≤ 8, i.e. `-ndiag 4` (or 1).
With only 2 k-points, `-npool 2` is the max useful pool count (48 ranks/pool);
the remaining parallelism must come from plane waves / bands / OpenMP.

**2. Out of memory.** **Decrease** `npool`. Pools duplicate the G-vector data, so
fewer pools means more ranks share (distribute) the plane-wave arrays, lowering
per-rank memory. (More nodes also helps.)

**3. Walltime hit.** Resubmit with **`recover = .true.`** in `&inputph`; `ph.x`
resumes from the last finished q/representation. You must **not** delete the
`.save`, `_ph0`, or `save` directories in the job script — that is the recover
state.

**4. The right machine.** 2-atom SCF → workstation, a few cores (k-point pools).
200-atom relaxation → HPC, plane-wave + band parallel (+ OpenMP), few k-points.
Phonon dispersion → HPC, `-npool` over k-points within each q, and `ph.x` images
over q/irreps if available.

</details>

---

## Congratulations 🎉

You've finished the tutorial. From a first SCF on GaAs to its band structure,
density of states, charge density, magnetism, and full phonon dispersion — with
the convergence discipline and HPC workflow to do it for real. The
[introduction](../intro) has the full chapter map; the runnable inputs live in
[`code/`](https://github.com/chaewoon11/qe-tutorial/tree/master/code). Go compute
something.
