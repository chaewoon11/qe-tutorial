#!/bin/bash
#PBS -N nbse2_ildos2
#PBS -V
#PBS -q debug
#PBS -A qe
#PBS -l select=2:ncpus=68:mpiprocs=68:ompthreads=1
#PBS -l walltime=02:00:00
#PBS -o ildos2.log
#PBS -e ildos2.err
cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh
module purge
module load intel/oneapi_21.2 impi/oneapi_21.2 craype-mic-knl
QE=$HOME/qe-7.5-knl/q-e-qe-7.5/bin
NP=$(wc -l < $PBS_NODEFILE)
echo "nscf(smearing) $(date)"
mpirun -np $NP $QE/pw.x -npool 4 < nbse2.3x3.cdw.nscf_sm.in > nbse2.3x3.cdw.nscf_sm.out
grep -q "JOB DONE" nbse2.3x3.cdw.nscf_sm.out && echo "nscf ok" || { echo FAIL; tail -15 nbse2.3x3.cdw.nscf_sm.out; exit 1; }
echo "Fermi: $(grep -i fermi nbse2.3x3.cdw.nscf_sm.out | tail -1)"
echo "pp ILDOS $(date)"
mpirun -np $NP $QE/pp.x < nbse2.3x3.cdw.ildos.in > nbse2.3x3.cdw.ildos.out
grep -q "Plot Type\|JOB DONE" nbse2.3x3.cdw.ildos.out && echo "pp ok" || { echo FAIL; tail -20 nbse2.3x3.cdw.ildos.out; }
ls -la nbse2_cdw.ildos.cube 2>/dev/null
