#!/bin/bash
# we should be done within one hour
#PBS -l walltime=2:00:00
# the number of cores we need (use procs=N, or nodes=A:ppn=B, with N=A*B)
#PBS -l procs=100
#PBS -q hpc
#PBS -e exec.err
#PBS -o exec.out
#PBS -M tntr@dtu.dk
#PBS -m abe

# change into work directory
cd $PBS_O_WORKDIR

size="100 200 500 800 1000 1500 2000 3000 4000"
thread="1 4 8 12 16 20 24 30 40 60 100"

module load mpi/gcc

OUT="out_${D}_${PBS_JOBNAME}.${PBS_JOBID}.out"


for t in $thread
do
	for s in $size
	do
		mpirun -n $t image_seg -file Data/gray.bmp -iters 400 -size $s >> "LOG/$OUT"
	done
done
