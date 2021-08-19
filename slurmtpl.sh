#!/bin/bash
#SBATCH --job-name="<jobname>"
#SBATCH --partition=eng-instruction
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 02:00:00
#SBATCH --output=<fullpath>.out
#SBATCH --error=<fullpath>.err