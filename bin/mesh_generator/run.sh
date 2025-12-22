#!/bin/bash

#SBATCH --job-name=mesh_generator
#SBATCH --partition=debug
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --mail-type=end
#SBATCH --mail-user=youremail

module purge
module load gcc

g++ mesh_generator.cpp
./a.out
rm a.out