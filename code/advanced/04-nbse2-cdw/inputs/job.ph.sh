#!/bin/bash
#PBS -N nbse2_ph
#PBS -V
#PBS -q debug
#PBS -A qe
#PBS -l select=2:ncpus=68:mpiprocs=68:ompthreads=1
#PBS -l walltime=12:00:00
#PBS -o ph.log
#PBS -e ph.err

cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi_21.2 impi/oneapi_21.2 craype-mic-knl
QE=$HOME/qe-7.5-knl/q-e-qe-7.5/bin
NP=$(wc -l < $PBS_NODEFILE)

echo "=== NbSe2 phonon chain START $(date) NP=$NP ==="

echo "--- [1/4] scf ---"
mpirun -np $NP $QE/pw.x -npool 8 < nbse2.scf.in > nbse2.scf.out
grep -q "JOB DONE" nbse2.scf.out || { echo "SCF FAILED"; tail -30 nbse2.scf.out; exit 1; }

echo "--- [2/4] ph.x (ldisp 6x6x1) ---"
mpirun -np $NP $QE/ph.x -npool 8 < nbse2.ph.in > nbse2.ph.out
grep -q "JOB DONE" nbse2.ph.out || { echo "PH FAILED"; tail -40 nbse2.ph.out; exit 1; }

echo "--- [3/4] q2r ---"
mpirun -np 1 $QE/q2r.x < nbse2.q2r.in > nbse2.q2r.out

echo "--- [4/4] matdyn ---"
mpirun -np 1 $QE/matdyn.x < nbse2.matdyn.in > nbse2.matdyn.out

echo "=== chain END $(date) ==="
echo "--- gamma (q=1) frequencies ---"
grep -A12 "q = " nbse2.ph.out | head -30
