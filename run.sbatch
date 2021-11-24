#!/bin/bash
#SBATCH -J yuanhui
#SBATCH --comment=NCCouple
#SBATCH -n 10
#SBATCH -N 6
#SBATCH --ntasks-per-node=5
#SBATCH -p normal
#SBATCH -o %j
#SBATCH -e %j

mpirun -n 6 ./NCCouple
