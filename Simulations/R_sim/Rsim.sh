#!/bin/bash
 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=04:00:00
#SBATCH --job-name=Rsim_hd41
#SBATCH --output=Rsim_hd41.out
#SBATCH --mail-type=FAIL

module load r/4.5.0
export PATH=$HOME/local/bin:$PATH
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

cd $SCRATCH/copula2/R_sim/

Rscript Gumbel-hd41.R
