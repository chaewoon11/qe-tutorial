# Chapter 0 — Exercise solutions

### 1. Convergence behaviour (`conv_thr = 1.0d-4`)

- **3 iterations** (vs. 6 with `1.0d-8`) — a looser target is reached sooner.
- Total energy: **−15.76389664 Ry**, vs. **−15.76393891 Ry** at `1.0d-8`.
- Difference: **~4×10⁻⁵ Ry ≈ 0.6 meV**. Tiny here, but the lesson is that the
  total energy is only as converged as `conv_thr` lets it be. For *energy
  differences* you want both runs converged to well below the difference you
  care about — `1.0d-6` Ry or tighter is the usual safe default.

### 2. The mixing knob (`mixing_beta = 0.3`)

- **14 iterations** (vs. 6 at `beta = 0.7`) — converges to essentially the same
  energy (**−15.76393888 Ry**), just slower.
- `mixing_beta` controls how much of the *new* density is blended into the old
  one each cycle. Large beta = bigger steps = faster *when it works*, but can
  oscillate or diverge for "hard" systems (metals, magnetic materials, big
  cells). A **smaller beta takes more, smaller, safer steps** — the standard fix
  when an SCF won't converge. Silicon is easy, so 0.7 is fine and faster.

### 3. Break it on purpose

With a bad pseudopotential filename, `pw.x` stops almost immediately:

```
     Error in routine readpp (1):
     file .../Si_does_not_exist.upf not found

     stopping ...
```

QE prints errors in a `%%%%`-bordered block naming the **routine** (`readpp`)
and the reason, then `stopping ...`. There is **no "JOB DONE"**. First debugging
move on any failed run: read the last ~20 lines of the output.

### 4. Read the geometry

```bash
grep -E "unit-cell volume|number of Kohn-Sham states" my.scf.out
```
```
     unit-cell volume          =     270.0114 (a.u.)^3
     number of Kohn-Sham states=            4
```

Silicon contributes **4 valence electrons per atom × 2 atoms = 8 electrons**.
This is a (non-magnetic) insulator with doubly-occupied states, so it needs
8 / 2 = **4 occupied bands** — hence 4 KS states by default. (When we add empty
bands for a band structure in Chapter 5, we'll raise this with `nbnd`.)
