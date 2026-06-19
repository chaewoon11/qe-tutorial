# Chapter 7 — charge density & visualization

```bash
pw.x < gaas.scf.in    > scf.out       # self-consistent density
pp.x < gaas.pp_110.in > pp_110.out    # 2D (1-10) slice -> gaas.rho_110.dat
pp.x < gaas.pp_xsf.in > pp_xsf.out    # 3D grid         -> gaas.rho.xsf (XCrySDen/VESTA)
```

- `gaas.rho_110.dat` — 2D charge-density slice (x, y, rho; bohr), plotted by
  `../../notebooks/07-charge-density.ipynb`.
- `gaas.rho.xsf`     — 3D density for isosurface visualization.
