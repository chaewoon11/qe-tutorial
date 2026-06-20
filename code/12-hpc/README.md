# Chapter 12 — running QE at scale

- `gaas.scf.in`     — a light GaAs SCF (ecut 50/400, 6×6×6 k, conv_thr 1e-8) used
  for the parallel scaling test.
- `strong_scan.dat` — wall time vs MPI ranks (`-npool = ranks`, OMP_NUM_THREADS=1).
- `npool_scan.dat`  — wall time vs `-npool` at fixed 8 ranks.

Strong-scaling test (vesper, OpenMPI):
```bash
export OMP_NUM_THREADS=1
for n in 1 2 4 8 16; do mpirun -np $n pw.x -npool $n < gaas.scf.in > scf_$n.out; done
```
Plotted by `../../notebooks/12-scaling.ipynb`. The 2-atom cell saturates by ~8
cores (speedup ~2.9×) — small jobs don't scale far.
