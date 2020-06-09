#! /bin/bash

#RES=$(sbatch slurm_gentrain.py) && sbatch --dependency=afterok:${RES##* } slurm_genmodel.py
# first job - no dependencies
jid1=$(sbatch slurm_gentrain.py | cut -f 4 -d' ')
echo "jid1 $jid1"

jid2=$(sbatch  --dependency=afterok:$jid1 slurm_genmodel.py | cut -f 4 -d' ')
echo "jid2 $jid2"

jid3=$(sbatch  --dependency=afterok:$jid2 slurm_selectmdl.py | cut -f 4 -d' ')
echo "jid3 $jid3"
