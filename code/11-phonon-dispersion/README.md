# Chapter 11 — phonon dispersion (GaAs)

Run on HPC (Nurion debug queue). From the converged SCF (ch10 settings):

```bash
pw.x     < gaas.scf.in        > gaas.scf.out        # density
ph.x     < gaas.ph.disp.in    > gaas.ph.disp.out    # DFPT on 4x4x4 q-grid (heavy)
q2r.x    < gaas.q2r.in        > gaas.q2r.out        # -> gaas.fc (force constants + Z*, eps)
matdyn.x < gaas.matdyn.in     > gaas.matdyn.out     # -> gaas.freq.gp  (dispersion, L-G-X-W-K-G)
matdyn.x < gaas.matdyn.dos.in > gaas.matdyn.dos.out # -> gaas.phdos   (phonon DOS, 24^3)
```

`ref/` holds the verified outputs; `gaas.freq.gp` + `gaas.phdos` are plotted by
`../../notebooks/11-phonon-dispersion.ipynb`.

Results (4x4x4 q-grid): Γ TO 269.7 cm⁻¹, LO 284.0 cm⁻¹ (LO-TO split ~14 cm⁻¹),
acoustic → 0 at Γ (ASR via matdyn). PBE LO-TO split underestimated (exp ~24)
because ε∞ is over-large.
