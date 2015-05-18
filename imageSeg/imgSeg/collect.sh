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
#PBS -l procs=8
#
# change into work directory
cd $PBS_O_WORKDIR


# load the right compiler and MPI support module
#
module load studio
module load mpi/studio

PNAME=`echo $PROG | sed 's/ .*$//'` 
echo $PNAME
# STR=`basename $PNAME`
# echo $STR
# JID=`echo $PBS_JOBID | sed 's/\..*$//'`

OUT="out_${D}_${PBS_JOBNAME}.${PBS_JOBID}.out"

# start collect with $PROG, using OpenMPI (aka CT8.1)
#
collect -M OPENMPI -o $OUT.er mpirun -np 6 -- ./image_seg -size 1000 -iters 400
