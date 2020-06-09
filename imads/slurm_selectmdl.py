#!/bin/env python

#SBATCH --mem=5G  
#SBATCH --mail-type=END
#SBATCH --mail-user=vm76@duke.edu

import os, sys
sys.path.append(os.getcwd())
import pickle
import libsvm.svmutil as svmutil

import imads_train as mt
from inputdict import param


if __name__ == "__main__":
    train = pickle.load(open('%s/imadstrain_w%d.pickle' % (param["outdir"],param["width"]), "rb" ))
    paramfiles = [fn for fn in os.listdir(param["outdir"]) if fn.startswith("parm_w%d_"%param["width"])]
    allparams = {}
    for fn in paramfiles:
        par = pickle.load(open('%s/%s' % (param["outdir"],fn), "rb" ))
        for core in par:
            if core not in allparams:
                allparams[core] = []
            allparams[core].append(par[core])

    model_fname =  '%s_w%s' % (param["tfname"], param["width"])
    param_log = ""
    for core in allparams:
        best_dict = max(allparams[core], key = lambda p:p["avg_scc"])
        model = mt.generate_svm_model(train[core], best_dict["params"], param["kmers"])
        svmutil.svm_save_model('%s/%s_%s.model' % (param["outdir"],model_fname,core), model)
        param_log += "%s: %s\n" % (core,str(best_dict))
    with open("%s/%s.log" % (param["outdir"],model_fname), 'w') as f:
        f.write(param_log)
