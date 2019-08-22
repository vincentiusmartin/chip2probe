#!/bin/bash
#SBATCH --job-name download_chip 
#SBATCH --mail-user th184@duke.edu
#SBATCH --mail-type ALL
#SBATCH --time 12:00:00
#SBATCH -c 3
#SBATCH --output "downloadchip".out
#SBATCH --error "downloadchip".err 
#SBATCH --mem=5G
file=$1
python3 $file
