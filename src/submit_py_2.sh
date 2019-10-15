#!/bin/bash
#SBATCH --job-name call_peaks_2 
#SBATCH --mail-user th184@duke.edu
#SBATCH --mail-type ALL
#SBATCH -c 3
#SBATCH --output "callpeaks_2".out
#SBATCH --error "callpeaks_2".err 
#SBATCH --mem=5G
file=$1
python3 $file
