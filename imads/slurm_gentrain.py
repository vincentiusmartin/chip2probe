#!/bin/env python

#SBATCH --mem=10G

import sys,os
sys.path.append(os.getcwd())

import pandas as pd
import pickle

import imads_train as mt
from inputdict import param

# sbatch --dependency=afterok:854 second_job.slurm
if __name__ == "__main__":

    cores_centered = mt.init_train_matrix(param)
    pickle.dump(cores_centered, open('%s/imadstrain_w%d.pickle' % (param["outdir"],param["width"]), 'wb'))
