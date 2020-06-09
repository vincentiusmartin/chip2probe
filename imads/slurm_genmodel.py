#!/bin/env python

#SBATCH -o outdir/vm76_%A_%a.out
#SBATCH -e outdir/vm76_%A_%a.err

#SBATCH --array=0-503%127
#SBATCH --mem=5G

import os, sys
sys.path.append(os.getcwd())
import pickle
import itertools

import imads_train as mt
from inputdict import param

def test_param_comb(traindata, param_dict, numfold=10, kmers=[1,2,3]):
    param_log = ""
    run_params = {}
    for core in traindata:
        #print("Working for core: %s" % core)
        run_params[core] = mt.run_kfold(param_dict, rows=traindata[core], numfold=numfold, kmers=kmers)
    return run_params

if __name__ == "__main__":
    aid = int(os.environ['SLURM_ARRAY_TASK_ID'])
    print("Array id: %d" % aid)

    train = pickle.load(open('%s/imadstrain_w%d.pickle' % (param["outdir"],param["width"]), "rb" ))
    # Make list of params for grid search
    paramkeys = list(param["grid"].keys())
    combinations = itertools.product(*param["grid"].values())
    combinations = [{paramkeys[i]:c[i] for i in range(len(c))} for c in combinations]
    print(len(combinations))
    selected_comb = combinations[aid]
    print(selected_comb)

    parm = test_param_comb(train, selected_comb, param["numfold"], param["kmers"])
    pickle.dump(parm, open('%s/parm_w%d_%d.pickle' % (param["outdir"],param["width"],aid), 'wb'))