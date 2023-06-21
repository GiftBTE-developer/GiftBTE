#!/bin/bash

#SBATCH --job-name=openmpitest
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load openmpi/4.1.1-gcc-8.3.1

mpirun -np 16 ../../BTE_CPU
