#!/bin/bash
#SBATCH --job-name="protoubiq"
#SBATCH --partition=eng-instruction
#SBATCH --nodes=10
#SBATCH --export=ALL
#SBATCH --ntasks-per-node=20
#SBATCH -t 04:00:00
export SLURM_NODEFILE=`generate_pbs_nodefile`
./julia --machinefile $SLURM_NODEFILE protoubiq.jl -i "test/*.txt" -e "ls {{}}" -s ".ls"