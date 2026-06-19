# Chapter 2 — cutoff convergence data

Total energy of GaAs as the plane-wave cutoffs are varied (all other settings as
in `../00-first-scf/inputs/gaas.scf.in`, with `conv_thr = 1.0d-10` for speed):

- `ecutwfc_conv.dat` — vary `ecutwfc` (20–80 Ry), `ecutrho = 400` fixed.
- `ecutrho_conv.dat` — vary `ecutrho` (200–600 Ry), `ecutwfc = 50` fixed.

Each file: two columns, `cutoff(Ry)  total_energy(Ry)`. Plotted by
`../../notebooks/02-cutoff-convergence.ipynb` as convergence curves in meV/atom.

Regenerate by looping the SCF input over cutoffs and grepping the `!  total
energy` line (see the loop in the Chapter 2 page).
