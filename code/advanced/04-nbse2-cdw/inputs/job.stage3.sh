#!/bin/bash
#PBS -N nbse2_stage3
#PBS -V
#PBS -q debug
#PBS -A qe
#PBS -l select=2:ncpus=68:mpiprocs=68:ompthreads=1
#PBS -l walltime=08:00:00
#PBS -o stage3.log
#PBS -e stage3.err

cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi_21.2 impi/oneapi_21.2 craype-mic-knl
QE=$HOME/qe-7.5-knl/q-e-qe-7.5/bin
NP=$(wc -l < $PBS_NODEFILE)
run(){ echo "--- $1 $(date +%T) ---"; mpirun -np $NP $QE/$2 -npool 4 < $3 > $4; grep -q "JOB DONE" $4 && echo ok || { echo "FAIL $3"; tail -15 $4; }; }

echo "=== Stage 3 START $(date) NP=$NP ==="
run "CDW scf"   pw.x  nbse2.3x3.cdw.scf.in   nbse2.3x3.cdw.scf.out
run "CDW pp"    pp.x  nbse2.3x3.cdw.pp.in    nbse2.3x3.cdw.pp.out
run "CDW nscf"  pw.x  nbse2.3x3.cdw.nscf.in  nbse2.3x3.cdw.nscf.out
mpirun -np 1 $QE/dos.x < nbse2.3x3.cdw.dos.in > nbse2.3x3.cdw.dos.out
run "ref scf"   pw.x  nbse2.3x3.ref.scf.in   nbse2.3x3.ref.scf.out
run "ref nscf"  pw.x  nbse2.3x3.ref.nscf.in  nbse2.3x3.ref.nscf.out
mpirun -np 1 $QE/dos.x < nbse2.3x3.ref.dos.in > nbse2.3x3.ref.dos.out
echo "=== END $(date) ==="
echo "Fermi (cdw): $(grep -i fermi nbse2.3x3.cdw.nscf.out | tail -1)"
echo "Fermi (ref): $(grep -i fermi nbse2.3x3.ref.nscf.out | tail -1)"
ls -la nbse2_cdw.charge.cube nbse2_cdw.dos nbse2_ref.dos 2>/dev/null
