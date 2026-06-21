#!/usr/bin/env python3
"""Freeze the q_CDW=(1/3,0) soft mode of 1H-NbSe2 into a 3x3 supercell.

Reads the soft-mode eigenvector (polarization vector) that ph.x already printed
in the q-point dynamical-matrix file (ref/nbse2.dyn8), and builds the modulated
3x3 supercell as either a single-q (1Q, stripe) or triple-q (3Q, triangular)
distortion. The 3Q pattern is the experimental NbSe2 ground state.

Output: QE relax inputs for the undistorted and distorted 3x3 cells.
No external packages beyond numpy.
"""
import numpy as np

# ---- relaxed 1H-NbSe2 primitive cell (from vc-relax) ----
A0 = 3.472867181                       # in-plane lattice constant (Ang)
C  = 22.0                              # vacuum (Ang)
a1 = np.array([ A0,        0.0,      0.0])
a2 = np.array([-A0/2, A0*np.sqrt(3)/2, 0.0])
a3 = np.array([ 0.0,       0.0,      C ])
# fractional atomic positions in the primitive cell
zSe = 0.076405410
frac = np.array([[1/3, 2/3, 0.0],      # Nb
                 [2/3, 1/3,  zSe],     # Se
                 [2/3, 1/3, -zSe]])    # Se
species = ['Nb', 'Se', 'Se']
mass = {'Nb': 92.90638, 'Se': 78.971}

# soft-mode polarization vectors e_i (Cartesian, complex) from ref/nbse2.dyn8, freq(1)
# lines: (Re Im) for x,y,z per atom
e0 = np.array([
    [-0.646031-0.300493j, -0.372988-0.173994j,  0.0+0.0j      ],  # Nb
    [ 0.118763-0.097775j,  0.068532-0.056520j, -0.360446+0.0j ],  # Se
    [ 0.118763-0.097775j,  0.068532-0.056520j,  0.360446+0.0j ],  # Se
])
# q_CDW in crystal (reduced) coordinates
q0_cryst = np.array([1/3, 0.0, 0.0])

# ---- reciprocal lattice (2pi included) ----
B = 2*np.pi*np.linalg.inv(np.array([a1, a2, a3]).T).T   # rows b1,b2,b3
def cryst_to_cart_q(qc):
    return qc @ B                                       # Cartesian q (1/Ang), incl 2pi

def rot120(n):
    """Cartesian rotation by 120*n degrees about z."""
    th = np.deg2rad(120*n)
    c, s = np.cos(th), np.sin(th)
    R = np.eye(3); R[0,0]=c; R[0,1]=-s; R[1,0]=s; R[1,1]=c
    return R

def build(mode='3Q', amp=0.10):
    """Return list of (species, cart_pos) for the modulated 3x3 supercell.
    amp = target max displacement in Ang."""
    members = [0] if mode == '1Q' else [0, 1, 2]
    # rotated (q, e) for each star member
    qs = [cryst_to_cart_q(rot120(m) @ q0_cryst) for m in members]   # note: rotate q in cart
    es = [ (rot120(m) @ e0.T).T for m in members ]                  # rotate each atom's e
    # mass-weight: physical displacement ~ e / sqrt(M)
    inv_sqrt_m = np.array([1/np.sqrt(mass[s]) for s in species])[:,None]

    raw = []   # (species, base_cart, disp_unscaled)
    for n1 in range(3):
        for n2 in range(3):
            T = n1*a1 + n2*a2                       # cell origin (Cartesian)
            for i,(sp,fp) in enumerate(zip(species, frac)):
                base = fp @ np.array([a1,a2,a3])    # atom Cartesian in primitive
                pos  = base + T
                u = np.zeros(3, complex)
                for q,e in zip(qs, es):
                    u += e[i] * np.exp(1j * np.dot(q, T))
                u = np.real(u) * inv_sqrt_m[i,0]
                raw.append([sp, pos, u])
    # scale so the largest displacement equals amp
    umax = max(np.linalg.norm(r[2]) for r in raw)
    scale = amp / umax
    atoms = [[r[0], r[1] + scale*r[2]] for r in raw]
    undist = [[r[0], r[1]] for r in raw]
    return atoms, undist

# 3x3 supercell lattice
A1, A2, A3 = 3*a1, 3*a2, a3

HEADER = """&control
    calculation = '{calc}'
    prefix      = 'nbse2_3x3'
    outdir      = './out'
    pseudo_dir  = '../pseudos'
    verbosity   = 'high'
    forc_conv_thr = 1.0d-5
    etot_conv_thr = 1.0d-6
/
&system
    ibrav = 0
    nat = 27
    ntyp = 2
    ecutwfc = 60.0
    ecutrho = 480.0
    assume_isolated = '2D'
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.01
/
&electrons
    conv_thr = 1.0d-12
    mixing_beta = 0.7
/
&ions
    ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
  Nb 92.90638 Nb.pbe-spn-kjpaw_psl.1.0.0.UPF
  Se 78.971   Se.pbe-n-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
{cell}
ATOMIC_POSITIONS angstrom
{pos}K_POINTS automatic
  4 4 1 0 0 0
"""

def write_input(fname, atoms, calc):
    cell = "\n".join("  %15.9f %15.9f %15.9f" % tuple(v) for v in (A1,A2,A3))
    pos  = "".join("  %-3s %15.9f %15.9f %15.9f\n" % (s, *p) for s,p in atoms)
    with open(fname, 'w') as f:
        f.write(HEADER.format(calc=calc, cell=cell, pos=pos))
    print("wrote", fname, "(%d atoms)" % len(atoms))

if __name__ == '__main__':
    atoms3Q, undist = build('3Q', amp=0.10)
    atoms1Q, _      = build('1Q', amp=0.10)
    write_input('nbse2.3x3.scf.in',     undist,  'scf')     # undistorted reference
    write_input('nbse2.3x3.relax.in',   atoms3Q, 'relax')   # 3Q triangular freeze-in
    write_input('nbse2.3x3.relax1Q.in', atoms1Q, 'relax')   # 1Q stripe freeze-in
    # report displacement character
    du = np.array([a[1]-u[1] for a,u in zip(atoms3Q, undist)])
    nb = du[0::3]; se = du[1::3].tolist()+du[2::3].tolist()
    print("3Q max |u_Nb| = %.3f A, max |u_Se| = %.3f A" %
          (np.linalg.norm(nb,axis=1).max(), max(np.linalg.norm(s) for s in se)))
