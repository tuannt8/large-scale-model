#!/bin/sh
#PBS -N week2
#PBS -o $PBS_JOBNAME.$PBS_JOBID.out
#PBS -e $PBS_JOBNAME.$PBS_JOBID.err
#PBS -l nodes=4:ppn=4
#PBS -q hpc

TARGET="gol"
NB_PROCS=10

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

module load mpi
MPIRUN=mpirun
$MPIRUN -n $NB_PROCS gol_block 10000	

