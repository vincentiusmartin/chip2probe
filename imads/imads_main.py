import pandas as pd
import itertools

from inputdict import param
import imads_train as mt

import libsvm.svmutil as svmutil


if __name__ == "__main__":
    data = pd.read_csv(param["pbmdata"], sep="\t", index_col="ID_REF")
    # just take the bound column, ignore the negative control
    df = data.loc[data["ID"] == "Bound"][[param["column_train"],"Sequence"]]

    # 1. Generate training data

    cores_centered = mt.gen_seqwcore(df.values.tolist(), param["width"], param["corelist"], corepos=param["corepos"])

    if param["logit"]:
        cores_cent = {k: [(mt.logit_score(val),seq) for (val, seq) in cores_centered[k]] for k in cores_centered}
    else:
        cores_cent = cores_centered

    # 2. Generate models

    paramkeys = list(param["grid"].keys())
    combinations = itertools.product(*param["grid"].values())
    combinations = [{paramkeys[i]:c[i] for i in range(len(c))} for c in combinations]

    paramscores = []
    for comb in combinations:
        parm = mt.test_param_comb(cores_cent, comb, param["numfold"], param["kmers"])
        paramscores.append(parm)

    # 3. Select best model

    allparams = {}
    for par in paramscores:
        for core in par:
            if core not in allparams:
                allparams[core] = []
            allparams[core].append(par[core])

    model_fname =  '%s_w%s' % (param["tfname"], param["width"])
    param_log = ""
    outdir = "%s/"%param["outdir"] if param["outdir"]  else ""
    for core in allparams:
        best_dict = max(allparams[core], key = lambda p:p["avg_scc"])
        model = mt.generate_svm_model(cores_cent[core], best_dict["params"], param["kmers"])
        svmutil.svm_save_model('%s%s_%s.model' % (outdir,model_fname,core), model)
        param_log += "%s: %s\n" % (core,str(best_dict))
    with open("%s%s.log" % (outdir,model_fname), 'w') as f:
        f.write(param_log)
