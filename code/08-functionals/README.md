# Chapter 8 — XC functionals & DFT+U

## GaAs LDA vs PBE (equation of state)
- `gaas.lda.scf.in` — GaAs scf with LDA (PZ) pseudopotentials; loop celldm(1) for the EOS.
- `eos_lda.dat`, `eos_pbe.dat` — total energy vs celldm(1) for each functional.
- PBE data reuses Chapter 4. Plotted by `../../notebooks/08-lda-vs-pbe.ipynb`
  (LDA a0=5.612 Å, PBE a0=5.756 Å, exp 5.653 Å).

## NiO DFT+U (AFM type-II)
```bash
pw.x < nio.pbe.scf.in  > nio.pbe.out     # HUBBARD U = 0.0001 (≈ PBE)
pw.x < nio.pbeu.scf.in > nio.pbeu.out    # HUBBARD U = 5.0 eV
```
- 4-atom AFM cell (Ni1 up, Ni2 down, 2 O), nspin=2.
- PBE: gap 1.0 eV, moment ~1.5 μB/Ni;  PBE+U(5 eV): gap 2.65 eV, moment ~1.7 μB/Ni.
