# Chapter 9 — spin & magnetism (bcc Fe)

```bash
pw.x  < fe.fm.scf.in > fe.fm.out     # ferromagnetic (nspin=2) -> moment 2.23 μB
pw.x  < fe.nm.scf.in > fe.nm.out     # nonmagnetic (nspin=1); FM is 0.57 eV/atom lower
pw.x  < fe.nscf.in   > fe.nscf.out   # dense tetrahedra grid for DOS
dos.x < fe.dos.in    > fe.dos.out    # spin-resolved DOS -> fe.dos
```

Note: dos.x Emin/Emax are *absolute* energies; Fe's E_F ≈ 17.4 eV, so the window
is 8–25 eV. `fe.dos` columns: E, dos_up, dos_dn, int_dos. Plotted by
`../../notebooks/09-fe-spin-dos.ipynb`.
