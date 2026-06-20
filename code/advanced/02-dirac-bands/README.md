# A2 — graphene Dirac bands

On the A1 relaxed cell:
```bash
pw.x    < graphene.scf.in      > scf.out        # density (conv_thr 1e-20)
pw.x    < graphene.bands.in    > bands.out      # bands along Γ-M-K-Γ (nbnd=12, conv_thr 1e-16)
bands.x < graphene.bands.pp.in > bands.pp.out   # -> graphene.bands.dat.gnu
```
E_F = -4.225 eV. The π/π* bands touch at the Dirac point K (gap ~0). Plotted by
`../../../notebooks/A2-dirac-bands.ipynb`.
