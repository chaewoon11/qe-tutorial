#!/bin/bash
#PBS -N nbse2_relax
#PBS -V
#PBS -q debug
#PBS -A qe
#PBS -l select=1:ncpus=68:mpiprocs=68:ompthreads=1
#PBS -l walltime=04:00:00
#PBS -o relax.log
#PBS -e relax.err

cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi_21.2 impi/oneapi_21.2 craype-mic-knl
QE=$HOME/qe-7.5-knl/q-e-qe-7.5/bin
NP=$(wc -l < $PBS_NODEFILE)

echo "=== NbSe2 vc-relax START $(date) NP=$NP ==="
mpirun -np $NP $QE/pw.x -npool 4 < nbse2.vc-relax.in > nbse2.vc-relax.out
grep -q "JOB DONE" nbse2.vc-relax.out && echo "DONE" || { echo "FAILED"; tail -30 nbse2.vc-relax.out; }
echo "=== END $(date) ==="
