# Chapter 6 — DOS & PDOS

Workflow (run in this directory):

```bash
pw.x      < gaas.scf.in   > scf.out      # density (8×8×8)
pw.x      < gaas.nscf.in  > nscf.out     # dense grid, tetrahedra, nbnd=24 (16×16×16)
dos.x     < gaas.dos.in   > dos.out      # -> gaas.dos (total DOS)
projwfc.x < gaas.proj.in  > proj.out     # -> gaas.pdos_* (projected DOS) + Lowdin charges
```

Committed outputs: `gaas.dos`, `gaas.pdos_atm#*(El)_wfc#*(l)`, `gaas.pdos_tot`.
Plotted by `../../notebooks/06-dos.ipynb` (total DOS + orbital-projected PDOS).
