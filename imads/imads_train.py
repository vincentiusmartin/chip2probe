# make svr model
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
import itertools
import ast
import sys
import libsvm.svmutil as svmutil
import os
import concurrent.futures as cc
import functools
import math
from tqdm import tqdm

#Q: why do we need to separate core?

def sparse_to_dense(sparse_matrix, totalfeat):
    dense = [float(sparse_matrix[i+1]) if i + 1 in sparse_matrix else 0.0 for i in range(totalfeat)]
    return dense

def count_dense_feat(lenseq, kmers):
    return sum([(lenseq+1-k) * pow(4,k) for k in kmers])

def libsvm_generate_matrix(seqlist, kmers, dense=False):
    """Generates the sparse matrix file from a list of sequences and their scores"""
    kmer_dict = {}
    for k in kmers:
        lst = ["".join(n) for n in itertools.product('ACGT', repeat=k)]
        kmer_dict[k] = {k: v + 1 for v, k in enumerate(lst)}
    scores = []
    features = []
    total_dense_feat = count_dense_feat(len(seqlist[0][1]), [1,2,3])
    for line in seqlist:
        score, seq = line
        scores.append(score)
        feat_pos = 0
        feat_dict = {}
        for k in kmers:
            hop = len(kmer_dict[k])
            for i in range(len(seq) - k + 1):
                kmer_seq = seq[i:i+k]
                feat_idx = feat_pos + kmer_dict[k][kmer_seq]
                feat_dict[feat_idx] = 1
                feat_pos += hop
        if dense:
            features.append(sparse_to_dense(feat_dict,total_dense_feat))
        else:
            features.append(feat_dict)
    return np.array(scores), np.array(features) # y, x

def gen_seqwcore(seqintensities, width, corelist, corepos="center"):
    """
    Generate core centered sequence and take mediann intensity if there are multiple matches
    :param seqintensities: list of list with [intensity, sequence]
    :param width:
    :param corelist:
    :param corepos: left, right, center
    :return: list of (score,sequence) with the sequence is of length 'width' with the core center
    """
    corelen = len(corelist[0])
    if not all(len(core) == corelen for core in corelist):
        raise ValueError('not all cores have same length!')
    core_dict = {}
    seqlen = len(seqintensities[0][1])
    if corepos == "left":
        s1 = int(math.ceil(0.5 * seqlen) - 0.5 * corelen)
        c1 = s1
    elif corepos == "right":
        s1 = int(math.ceil(0.5 * seqlen) - width + 0.5 * corelen)
        c1 = s1 + width - corelen
    else: #center
        s1 = int(math.ceil(0.5 * seqlen) - 0.5 * width)
        c1 = int(0.5 * seqlen - 0.5 * corelen)
    spos = (s1, s1 + width)
    cpos = (c1, c1 + corelen)
    # process each core separately and make sure that the list is unique
    for core in set(corelist):
        seq_wcore = [(score, seq[spos[0]:spos[1]]) for score, seq in seqintensities if seq[cpos[0]:cpos[1]] == core]
        seq_core_df = pd.DataFrame(seq_wcore, columns = ['score', 'seq'])
        agg = seq_core_df.groupby(['seq'], as_index=False).median()
        core_dict[core] = agg[["score","seq"]].values.tolist()
    return core_dict

def run_kfold(param_dict, rows, numfold, kmers=[1,2,3]):
    kf = KFold(numfold, shuffle=True)
    splitted = kf.split(rows)
    param_str = "-s 3 -b 1 -q " # epsilon-SVR, prob estimate true, quiet mode
    param_str += " ".join(["-{} {}".format(k,v) for k,v in param_dict.items()])
    params = svmutil.svm_parameter(param_str)
    scc_list = [] # r squared
    mse_list = []
    for train_idx, test_idx in splitted:
        train_list = [rows[i] for i in train_idx]
        test_list = [rows[i] for i in test_idx]

        y_train, x_train = libsvm_generate_matrix(train_list, kmers)
        y_test, x_test = libsvm_generate_matrix(test_list, kmers)

        train_prob  = svmutil.svm_problem(y_train, x_train)

        model = svmutil.svm_train(train_prob, params)
        #svmutil.svm_save_model('model_name.model', m)
        # y is only needed when we need the model performance
        p_label, p_acc, p_val = svmutil.svm_predict(y_test, x_test, model, options="-q")
        acc, mse, scc = p_acc
        scc_list.append(scc)
        mse_list.append(mse)
    avg_scc = sum(scc_list)/numfold
    avg_mse = sum(mse_list)/numfold
    return {"params":param_dict, "avg_scc": avg_scc, "avg_mse":avg_mse}

def test_param_comb(traindata, param_dict, numfold=10, kmers=[1,2,3]):
    param_log = ""
    run_params = {}
    for core in traindata:
        #print("Working for core: %s" % core)
        run_params[core] = run_kfold(param_dict, rows=traindata[core], numfold=numfold, kmers=kmers)
    return run_params

def generate_svm_model(rows, param, kmers=[1,2,3]):
    y, x = libsvm_generate_matrix(rows, kmers, dense=True)
    prob  = svmutil.svm_problem(y, x)
    # s = 3 -- epsilon-SVR
    param_str = "-s 3 -b 1 -q " # epsilon-SVR, prob estimate true, quiet mode
    param_str += " ".join(["-{} {}".format(k,v) for k,v in param.items()])
    params = svmutil.svm_parameter(param_str)
    model = svmutil.svm_train(prob, params)
    return model

def logit_score(p):
    # f(x) = 1 / ( 1 + exp(-x) )  to obtain only values between 0 and 1.
    return  np.log(p/(1.0-p))

def genmodel_gridsearch(cores_centered, param_dict, numfold=10, kmers=[1,2,3],
                        numworkers=os.cpu_count(), logit=True, tfname="",
                        modelwidth=20, outdir="", suffix=""):
    if logit:
        cores_cent = {k: [(logit_score(val),seq) for (val, seq) in cores_centered[k]] for k in cores_centered}
    else:
        cores_cent = cores_centered
    # Make list of params for grid search
    params = list(param_dict.keys())
    combinations = itertools.product(*param_dict.values())
    combinations = [{params[i]:c[i] for i in range(len(c))} for c in combinations]

    # ----- RUNNING CROSS VALIDATION -----
    param_log = ""
    model_fname =  '%s_w%s' % (tfname, modelwidth)
    if suffix:
        model_fname = "%s_%s" % (model_fname, suffix)
    outdir = "%s/"%outdir if outdir  else ""
    for core in cores_cent:
        print("Working for core: %s" % core)
        run_kfold_partial = functools.partial(run_kfold, rows=cores_cent[core], numfold=numfold, kmers=kmers)
        with cc.ProcessPoolExecutor(max_workers = numworkers) as executor:
            # TODO: update input to combinations to dictionary
            run_params = list(tqdm(executor.map(run_kfold_partial, combinations), total=len(combinations)))
        # Save the best model
        best_dict = max(run_params, key = lambda p:p["avg_scc"])

        model = generate_svm_model(cores_cent[core], best_dict["params"], kmers)
        svmutil.svm_save_model('%s%s_%s.model' % (outdir,model_fname,core), model)
        param_log += "%s: %s\n" % (core,str(best_dict))

    with open("%s%s.log" % (outdir,model_fname), 'w') as f:
        f.write(param_log)


def get_weight(libsvm_model, width, kmers=[1,2,3]):
    lendense = count_dense_feat(width, kmers)
    dense_vectors = np.array([sparse_to_dense(v, lendense) for v in mdl.get_SV()])
    coefs = np.transpose(np.array(mdl.get_sv_coef()))
    w = np.dot(coefs, dense_vectors)
    return w[0]

def get_feature_imp(weights, width):
    nucleotide = ["A", "C", "G", "T"]
    comb_dict = {}
    for k in kmers:
        comb_dict[k] = ["".join(n) for n in itertools.product('ACGT', repeat=k)]
    features = []
    for k in comb_dict:
        for i in range(1, width - k + 2):
            for comb in comb_dict[k]:
                features.append([comb, i])
    df = pd.DataFrame(features, columns = ["feature", "position"])
    df['weight'] = weights
    sorted = df.iloc[(-df['weight'].abs()).argsort()]
    return sorted

def explain_imp(impdf):
    df = pd.DataFrame(impdf)[["position", "weight"]]
    df["weight"] = abs(df["weight"])
    bypos = df.groupby("position")[["weight"]].sum().sort_values("weight", ascending=False)
    print(bypos)


"""
if __name__ == "__main__":
    pbmdata = "/Users/vincentiusmartin/Research/chip2gcPBM/imads/data/Combined_ets1_100nM_elk1_100nM_50nM_gabpa_100nM_log_normalized.txt"
    # normlized version
    column_train = "Ets1_100nM"
    kmers = [1,2,3]
    width = 20
    corelist = ["GGAA", "GGAT"] #, "GGAT"]
    tfname = "Ets1"
    param_dict = { #SVR feature
        "c": [0.05, 0.5, 1], # cost
        "g": [0.0005, 0.005, 0.01], # gamma for epsilon SVR
        "p": [0.01, 0.02, 0.1], # epsilon, linear don't use: 0.01
        "t": [2] # kernel type: 0:linear, 1: polynomial, 2: radial, 3: sigmoid, 4: precomputed
    }
    num_workers = os.cpu_count()
    numfold = 5

    # First, read pbmdata and generate kmer with their score to cross_centered
    data = pd.read_csv(pbmdata, sep="\t", index_col="ID_REF")
    df = data.loc[data["ID"] == "Bound"][[column_train,"Sequence"]]

    #for cp in ["center", "left", "right"]:
    cp="center"
    cores_centered = gen_seqwcore(df.values.tolist(), width, corelist, corepos=cp)

    mdl = svmutil.svm_load_model("model/w12/Ets1_w12_center_GGAA.model")
    # # mdl = generate_svm_model(cores_centered["GGAA"], {'c': 0.01, 'p': 0.1, 't': 0})
    w = get_weight(mdl, width, kmers)
    imp = get_feature_imp(w, width)
    explain_imp(imp)
    imp.to_csv("Ets1_w12_center_GGAA.model_imp.csv",index=False)

    #genmodel_gridsearch(cores_centered, param_dict, numfold, kmers, num_workers, logit=True, suffix=cp)
"""
