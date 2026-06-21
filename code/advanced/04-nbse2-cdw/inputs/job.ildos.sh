#!/bin/bash
#PBS -N nbse2_ildos
#PBS -V
#PBS -q debug
#PBS -A qe
#PBS -l select=1:ncpus=68:mpiprocs=68:ompthreads=1
#PBS -l walltime=01:00:00
#PBS -o ildos.log
#PBS -e ildos.err
cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi_21.2 impi/oneapi_21.2 craype-mic-knl
QE=$HOME/qe-7.5-knl/q-e-qe-7.5/bin
echo "ILDOS start $(date)"
mpirun -np 68 $QE/pp.x < nbse2.3x3.cdw.ildos.in > nbse2.3x3.cdw.ildos.out
grep -q "JOB DONE\|Plot Type" nbse2.3x3.cdw.ildos.out && echo ok || { echo FAIL; tail -20 nbse2.3x3.cdw.ildos.out; }
ls -la nbse2_cdw.ildos.cube 2>/dev/null
