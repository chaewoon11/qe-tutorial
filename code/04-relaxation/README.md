# Chapter 4 — structural optimization

- `gaas.scf.in`       — GaAs scf (8×8×8 k-grid); celldm(1) is edited in the EOS sweep.
- `gaas.vc-relax.in`  — variable-cell relaxation of GaAs (forces + stress on).
- `eos.dat`           — total energy vs celldm(1) (bohr), for the equation-of-state fit.

Plotted/fitted by `../../notebooks/04-equation-of-state.ipynb` (Birch–Murnaghan
→ a0 ≈ 5.756 Å, B0 ≈ 60 GPa). vc-relax converges to a0 ≈ 5.748 Å.
