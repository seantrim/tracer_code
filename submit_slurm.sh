#!/bin/bash

#SBATCH --account=def-samb71
#SBATCH -t 0-12:00
#SBATCH --mem=10000
#SBATCH -c 8
#SBATCH --output=std.out
#SBATCH --error=std.err
module load intel
ulimit -s unlimited

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_STACKSIZE=15m
echo "Job ID="$SLURM_JOB_ID
lscpu
./Makefile
./a.out
