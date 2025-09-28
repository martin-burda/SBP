#!/bin/bash
 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=00:30:00
#SBATCH --job-name=mixture-hd41
#SBATCH --output=mixture-hd41.out
#SBATCH --mail-type=FAIL

module load StdEnv/2023 intel/2023.2.1

cd $SCRATCH/copula2/fortran_sim/
ifort random.f90 sim_SBP.f90 -o sim_SBP.exe -O3

./sim_SBP.exe