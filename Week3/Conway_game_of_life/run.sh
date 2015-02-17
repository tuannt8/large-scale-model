#!/bin/sh
#PBS -N week2
#PBS -o $PBS_JOBNAME.$PBS_JOBID.out
#PBS -e $PBS_JOBNAME.$PBS_JOBID.err
#PBS -l nodes=4:ppn=4
#PBS -q hpc

TARGET="mat_vec_mul_nonblock mat_vec_mul_nonblock_row"

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE

module load mpi
MPIRUN=mpirun
$MPIRUN -n 10 mat_vec_mul_nonblock 70000
$MPIRUN -n 10  mat_vec_mul_nonblock_row 70000 

