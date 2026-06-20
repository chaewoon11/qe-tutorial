# A3 — bilayer hBN: vdW, bands & 2D phonons (LO-TO)

AA' (AAp) bilayer hBN, 4 atoms, vdW-DF (DFT-D3), assume_isolated='2D', ecutwfc=98.
HPC-grade reference data (full 12×12×1 DFPT). Workflow:

```bash
# bands:
pw.x    < bn.scf.in      > scf.out
pw.x    < bn.bands.in    > bands.out      # -> bn.bands.gnu (Γ-M-K-Γ)
bands.x < bn.bands.pp.in > bands.pp.out
# phonons (HPC):
pw.x     < bn.scf.in     > scf.out
ph.x     < bn.ph.in      > ph.out         # DFPT 12×12×1 q-grid, epsil (Born charges)
q2r.x    < bn.q2r.in     > q2r.out        # loto_2d=.true. -> bn.fc
matdyn.x < bn.matdyn.in  > matdyn.out     # loto_2d=.true. -> bn.freq.gp
```

Committed data: `bn.bands.gnu` (PBE indirect gap 4.53 eV, VBM~K / CBM~M) and
`bn.freq.gp` (phonon dispersion; E2g ~1352 cm⁻¹ with 2D LO-TO). Plotted by
`../../../notebooks/A3-bilayer-hbn.ipynb`.

Pseudopotentials: PSlibrary B/N PAW (PBE). The bundled reference data was computed
with equivalent B/N PAW potentials.
