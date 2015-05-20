#!/bin/bash
#
# collect.sh - a PBS script to start a collect session (Studio # analyzer)
#              for MPI programs (multi-node)
#
# submit as:  PROG=./your_prog qsub [options] collect.sh
#
# written by: Bernd Dammann <support@hpc.dtu.dk> , 03/2014
#
# name of the job
#PBS -N collect
# we need to forward the name of our program to debug
#PBS -v PROG
# keep the log files local
#PBS  -k oe
# we should be done within one hour
#PBS -l walltime=0:05:00
# the number of cores we need (use procs=N, or nodes=A:ppn=B, with N=A*B)
#PBS -l procs=64
#
# change into work directory
cd $PBS_O_WORKDIR


# load the right compiler and MPI support module
#
module load studio
module load mpi/studio

# OUT="out_${D}_${PBS_JOBNAME}.${PBS_JOBID}.out"
nproc="64"
SIZE="4000"
ITERS="100"
OUT="${nproc}_${SIZE}_${ITERS}_${PBS_JOBID}"

# start collect with $PROG, using OpenMPI (aka CT8.1)
#
collect -M OPENMPI -o $OUT.er mpirun -np $nproc -- ./image_seg_no_topo -size 4000 -iters $ITERS
