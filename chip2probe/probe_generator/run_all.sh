#!/bin/bash

#SBATCH -o vm_%A_%a.out
#SBATCH -e vm_%A_%a.err

#SBATCH --array=1%1
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=vm76@duke.edu

# TO RUN: sbatch run_all.sh

index=2 #${input[SLURM_ARRAY_TASK_ID]}
outdir="/Users/vincentiusmartin/Research/chip2gcPBM/result" #"/Users/vincentiusmartin/Research/chip2gcPBM/result" "/data/gordanlab/vincentius/cooperative_probe/result"

echo "python run_all.py $index -o $outdir"
python3 run_all.py $index -o $outdir
