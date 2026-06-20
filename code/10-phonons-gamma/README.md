# Chapter 10 — DFPT Γ phonons (GaAs)

Run (best on HPC — these are the first expensive jobs; we used Nurion debug queue):

```bash
pw.x < gaas.scf.in      > gaas.scf.out      # insulator SCF (ecut 70/560, 12^3 k, conv_thr 1e-20)
ph.x < gaas.ph.gamma.in > gaas.ph.out       # Γ phonons; epsil -> dielectric + Born charges; nmix_ph=12
```

`ref/` holds the verified outputs (gaas.scf.out, gaas.ph.out, gaas.dyn).

Converged results: acoustic at Γ −3.5 cm⁻¹ (raw, pre-ASR), TO 269.7 cm⁻¹,
ε∞ 17.8, Z* ±2.16. The acoustic frequency is the convergence gauge — the loose
ecut 50/400, 8^3 run gave 13.7 cm⁻¹ (under-converged); SCF conv_thr alone did not
help, the k-grid and cutoffs did.
