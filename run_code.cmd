#!/bin/bash

#SBATCH -t 100:00:00
#SBATCH --job-name=run
#SBATCH -p normal
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --output=run.txt
#SBATCH --ntasks-per-node=80

module load R

Rscript ./code/model_simulate.R