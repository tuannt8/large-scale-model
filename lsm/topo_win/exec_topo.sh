#!/bin/bash
# we should be done within one hour
#PBS -l walltime=01:00:00
# the number of cores we need (use procs=N, or nodes=A:ppn=B, with N=A*B)
#PBS -l procs=64
#PBS -q hpc
#PBS -e exec.err
#PBS -o exec.out
#PBS -M tntr@dtu.dk
#PBS -m abe

# change into work directory
cd $PBS_O_WORKDIR

module load mpi/gcc

OUT="out_${D}_${PBS_JOBNAME}.${PBS_JOBID}.out"

size="4000"
thread="63"
ITER="50000"
file="Data/gray.bmp"

for t in $thread
do
	for s in $size
	do
		time mpirun -n $t image_seg_topo -file $file -iters $ITER >> "LOG/$OUT"
		time mpirun -n $t image_seg_topo_win -file $file  -iters $ITER >> "LOG/$OUT"
		time mpirun -n $t image_seg_no_topo -file $file  -iters $ITER >> "LOG/$OUT"
	done
done
