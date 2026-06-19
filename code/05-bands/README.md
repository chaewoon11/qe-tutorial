# Chapter 5 — band structure

Three-step workflow (run in this directory):

```bash
pw.x    < gaas.scf.in       > scf.out        # self-consistent density
pw.x    < gaas.bands.in     > bands.out      # eigenvalues along L-G-X-W-K-G (nbnd=16)
bands.x < gaas.bands.pp.in  > bands.pp.out   # -> gaas.bands.dat.gnu
```

- `gaas.bands.dat.gnu` — committed band data (k-distance, energy), one band per block.
- `highsym.dat`        — high-symmetry tick positions + labels for the plot.

Plotted by `../../notebooks/05-band-structure.ipynb` (direct gap 0.51 eV at Γ).
