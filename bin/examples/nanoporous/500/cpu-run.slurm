#!/bin/bash

#SBATCH --job-name=openmpitest
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=4
#SBATCH -n 4
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module purge
module load openmpi/4.1.1-gcc-8.3.1

mpirun -np 4 ../../../BTE_CPU
