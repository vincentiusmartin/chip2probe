#!/bin/bash
#SBATCH --job-name call_peaks 
#SBATCH --mail-user th184@duke.edu
#SBATCH --mail-type ALL
#SBATCH -c 3
#SBATCH --output "callpeaks".out
#SBATCH --error "callpeaks".err 
#SBATCH --mem=5G
file=$1
python3 $file
