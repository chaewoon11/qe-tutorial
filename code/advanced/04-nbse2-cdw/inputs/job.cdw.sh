#!/bin/bash
#PBS -N nbse2_cdw3x3
#PBS -V
#PBS -q debug
#PBS -A qe
#PBS -l select=2:ncpus=68:mpiprocs=68:ompthreads=1
#PBS -l walltime=12:00:00
#PBS -o cdw.log
#PBS -e cdw.err

cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi_21.2 impi/oneapi_21.2 craype-mic-knl
QE=$HOME/qe-7.5-knl/q-e-qe-7.5/bin
NP=$(wc -l < $PBS_NODEFILE)

echo "=== NbSe2 3x3 CDW START $(date) NP=$NP ==="

echo "--- [1/3] undistorted 3x3 reference scf ---"
mpirun -np $NP $QE/pw.x -npool 4 < nbse2.3x3.scf.in > nbse2.3x3.scf.out
grep -q "JOB DONE" nbse2.3x3.scf.out || { echo "REF SCF FAILED"; tail -20 nbse2.3x3.scf.out; }

echo "--- [2/3] 3Q (triangular) freeze-in relax ---"
mpirun -np $NP $QE/pw.x -npool 4 < nbse2.3x3.relax.in > nbse2.3x3.relax.out
grep -q "JOB DONE" nbse2.3x3.relax.out || { echo "3Q RELAX FAILED"; tail -20 nbse2.3x3.relax.out; }

echo "--- [3/3] 1Q (stripe) freeze-in relax ---"
mpirun -np $NP $QE/pw.x -npool 4 < nbse2.3x3.relax1Q.in > nbse2.3x3.relax1Q.out
grep -q "JOB DONE" nbse2.3x3.relax1Q.out || { echo "1Q RELAX FAILED"; tail -20 nbse2.3x3.relax1Q.out; }

echo "=== END $(date) ==="
echo "--- energies (Ry) ---"
echo "ref(undist): $(grep '^!' nbse2.3x3.scf.out | tail -1)"
echo "3Q  final  : $(grep '^!' nbse2.3x3.relax.out | tail -1)"
echo "1Q  final  : $(grep '^!' nbse2.3x3.relax1Q.out | tail -1)"
