# A1 — building a 2D cell (graphene)

Professional workflow: converge (vacuum) -> relax -> reuse the relaxed cell.
conv_thr convention (2D track): scf 1e-20, relax 1e-18, nscf 1e-16, ph tr2_ph 1e-20.

```bash
# 1) vacuum convergence (k=12x12x1 fixed; loose conv for the test):
for c in 14 18 22 26 30; do
  sed "s/__C__/$c.0/" graphene.scf.tmpl > /tmp/g.in
  echo "$c $(grep '^!' <(pw.x < /tmp/g.in) | tail -1 | awk '{print $5}')"
done
# 2) relax the in-plane lattice (vc-relax, cell_dofree=2Dxy):
pw.x < graphene.relax.in > graphene.relax.out   # -> a0 = 2.4671 Ang
# 3) graphene.scf.in already holds the relaxed cell (conv_thr 1e-20) for reuse.
```

- `graphene.scf.in`   — converged + relaxed cell (a=2.4671 Å, c=20 Å, 2D, 12×12×1, mv smearing).
- `graphene.relax.in` — vc-relax (2Dxy) from experimental a.
- `graphene.scf.tmpl` — vacuum-sweep template (`__C__` placeholder).
- `vac.dat`           — vacuum convergence data; plotted by `../../../notebooks/A1-vacuum.ipynb`.

Gotcha: c ≲ 10 Å errors ("2D cutoff smaller than pseudo cutoff radius").
