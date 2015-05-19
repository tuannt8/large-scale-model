#!/bin/bash
# we should be done within one hour
#PBS -l walltime=00:1:00
# the number of cores we need (use procs=N, or nodes=A:ppn=B, with N=A*B)
#PBS -l procs=64
#PBS -q hpc
#PBS -e exec.err
#PBS -o exec.out
#PBS -M tntr@dtu.dk
# PBS -m abe

# change into work directory
cd $PBS_O_WORKDIR

size="4000"
thread="64"

module load mpi/gcc

OUT="out_${D}_${PBS_JOBNAME}.${PBS_JOBID}.out"


for t in $thread
do
	for s in $size
	do
		mpirun -n $t image_seg -size $s  -iters 200 >> "LOG/$OUT"
	done
done
