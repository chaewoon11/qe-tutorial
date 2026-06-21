#!/usr/bin/env python3
"""Generate Stage-3 inputs: charge density (pp.x) + DOS (nscf->dos.x) for the
undistorted and the relaxed 3x3 CDW cells of NbSe2."""
import numpy as np

a0 = 3.472867181
A1 = (3*a0, 0., 0.); A2 = (-3*a0/2, 3*a0*np.sqrt(3)/2, 0.); A3 = (0., 0., 22.)

def read_relax_final(fn):
    blk = open(fn).read().split('Begin final coordinates')[1].split('End final coordinates')[0]
    blk = blk.split('ATOMIC_POSITIONS')[1].splitlines()[1:]
    at = []
    for l in blk:
        t = l.split()
        if len(t) >= 4 and t[0] in ('Nb', 'Se'): at.append((t[0], float(t[1]), float(t[2]), float(t[3])))
    return at

def read_in_pos(fn):
    blk = open(fn).read().split('ATOMIC_POSITIONS angstrom')[1].split('K_POINTS')[0]
    at = []
    for l in blk.strip().splitlines():
        t = l.split()
        if t[0] in ('Nb', 'Se'): at.append((t[0], float(t[1]), float(t[2]), float(t[3])))
    return at

cdw = read_relax_final('ref/nbse2.3x3.relax.out')
ref = read_in_pos('nbse2.3x3.scf.in')
assert len(cdw) == 27 and len(ref) == 27

def cell_block():
    return ("CELL_PARAMETERS angstrom\n"
            "  %15.9f %15.9f %15.9f\n  %15.9f %15.9f %15.9f\n  %15.9f %15.9f %15.9f\n"
            % (*A1, *A2, *A3))

def pos_block(at):
    return "ATOMIC_POSITIONS angstrom\n" + "".join("  %-3s %15.9f %15.9f %15.9f\n" % a for a in at)

SPECIES = ("ATOMIC_SPECIES\n  Nb 92.90638 Nb.pbe-spn-kjpaw_psl.1.0.0.UPF\n"
           "  Se 78.971   Se.pbe-n-kjpaw_psl.1.0.0.UPF\n")

def pw(calc, prefix, outdir, at, kline, extra_sys=""):
    return f"""&control
    calculation = '{calc}'
    prefix      = '{prefix}'
    outdir      = '{outdir}'
    pseudo_dir  = '../../pseudos'
    verbosity   = 'high'
/
&system
    ibrav = 0
    nat = 27
    ntyp = 2
    ecutwfc = 60.0
    ecutrho = 480.0
    assume_isolated = '2D'
{extra_sys}/
&electrons
    conv_thr = 1.0d-10
    mixing_beta = 0.7
/
{SPECIES}{cell_block()}{pos_block(at)}{kline}"""

SM = "    occupations = 'smearing'\n    smearing = 'mv'\n    degauss = 0.01\n"
TT = "    occupations = 'tetrahedra'\n    nbnd = 250\n"

KSCF = "K_POINTS automatic\n  4 4 1 0 0 0\n"
KNSCF = "K_POINTS automatic\n  8 8 1 0 0 0\n"

# --- CDW (relaxed) ---
open('nbse2.3x3.cdw.scf.in','w').write(pw('scf','nbse2_cdw','./out_cdw',cdw,KSCF,SM))
open('nbse2.3x3.cdw.nscf.in','w').write(pw('nscf','nbse2_cdw','./out_cdw',cdw,KNSCF,TT))
# --- undistorted reference ---
open('nbse2.3x3.ref.scf.in','w').write(pw('scf','nbse2_ref','./out_ref',ref,KSCF,SM))
open('nbse2.3x3.ref.nscf.in','w').write(pw('nscf','nbse2_ref','./out_ref',ref,KNSCF,TT))

# --- dos.x ---
def dos(prefix, outdir, fil):
    return f"&dos\n  prefix = '{prefix}'\n  outdir = '{outdir}'\n  fildos = '{fil}'\n  emin = -8.0\n  emax = 6.0\n  deltae = 0.02\n/\n"
open('nbse2.3x3.cdw.dos.in','w').write(dos('nbse2_cdw','./out_cdw','nbse2_cdw.dos'))
open('nbse2.3x3.ref.dos.in','w').write(dos('nbse2_ref','./out_ref','nbse2_ref.dos'))

# --- pp.x: total valence charge of the CDW cell (3D cube) ---
open('nbse2.3x3.cdw.pp.in','w').write(
"""&inputpp
  prefix = 'nbse2_cdw'
  outdir = './out_cdw'
  plot_num = 0
  filplot = 'nbse2_cdw.chg'
/
&plot
  iflag = 3
  output_format = 6
  fileout = 'nbse2_cdw.charge.cube'
/
""")
print("wrote: cdw/ref scf+nscf+dos, cdw pp  (7 inputs)")
