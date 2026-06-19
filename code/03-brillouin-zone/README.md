# Chapter 3 — k-point convergence data

- `al.scf.in` — fcc aluminium, a one-atom metal with `mv` smearing.
- `gaas_kconv.dat` — GaAs total energy vs N (N×N×N grid), `occupations='fixed'`.
- `al_kconv.dat`   — Al total energy vs N, `occupations='smearing'`.

Both sweeps used `conv_thr=1.0d-10`. GaAs reuses `../00-first-scf/inputs/gaas.scf.in`
with the `K_POINTS` line edited. Plotted by `../../notebooks/03-kpoint-convergence.ipynb`.
