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

    # data = pd.read_csv(param["pbmdata"], sep="\t", index_col="ID_REF")
    #
    # # just take the bound column, ignore the negative control
    # bound_idxs = data[param["column_id"]].str.contains("Bound")
    # df = data[bound_idxs].reset_index()[[param["column_train"],"Sequence"]]
    #
    # cores_centered = mt.gen_seqwcore(df.values.tolist(), param["width"], param["corelist"], corepos=param["corepos"])
    #
    # if param["logit"]:
    #     cores_cent = {k: [(mt.logit_score(val),seq) for (val, seq) in cores_centered[k]] for k in cores_centered}
    # else:
    #     cores_cent = cores_centered
    cores_centered = mt.init_train_matrix(param)
    pickle.dump(cores_centered, open('%s/imadstrain_w%d.pickle' % (param["outdir"],param["width"]), 'wb'))
