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

# Wrapper function, accept param dictionary

def init_train_matrix(param):
    # make output directory if not exist
    if not os.path.exists(param["outdir"]):
        os.makedirs(param["outdir"])

    data = pd.read_csv(param["pbmdata"], sep="\t", index_col="ID_REF")
    df = pd.DataFrame(data[[param["column_id"],param["column_train"],"Sequence"]])
    if param["normalize"]:
        maxval, minval = df[param["column_train"]].max(), df[param["column_train"]].min()
        df[param["column_train"]] = df[param["column_train"]].apply(lambda x : (x - minval) / (maxval-minval))
    bound_idxs = df[param["column_id"]].str.contains("Bound")
    # just take the bound column, ignore the negative control
    df = df[bound_idxs].reset_index()[[param["column_train"],"Sequence"]]

    cores_centered = gen_seqwcore(df.values.tolist(), param["width"], param["corelist"], corepos=param["corepos"])

    # do the logistic transformation
    if param['logit']:
        cores_centered = {k: [(logit_score(val),seq) for (val, seq) in cores_centered[k]] for k in cores_centered}

    return cores_centered

def write_result(core_params, cores_centered, user_param):
    model_fname =  '%s_w%s' % (user_param['tfname'], user_param['width'])
    outdir = "%s/"%user_param['outdir'] if user_param['outdir']  else ""
    # if suffix:
    #     model_fname = "%s_%s" % (model_fname, suffix)
    param_log = ""
    pmlist = []
    for core in core_params:
        # Save the best model
        best_dict = max(core_params[core], key = lambda p:p["avg_scc"])

        # get predicted vs measured
        pm = predict_kfold(best_dict['params'], rows=cores_centered[core], numfold=user_param['numfold'], kmers=user_param['kmers'])
        pm = pd.DataFrame(pm)
        pm["core"] = core
        pmlist.append(pm)

        model = generate_svm_model(cores_centered[core], best_dict["params"], user_param['kmers'])
        svmutil.svm_save_model('%s%s_%s.model' % (outdir,model_fname,core), model)
        param_log += "%s: %s\n" % (core,str(best_dict))

    pm_df = pd.concat(pmlist)
    rsq_all = pm_df["measured"].corr(pm_df["predicted"])**2
    param_log += "R²: %s\n" % rsq_all

    with open("%s%s.log" % (outdir,model_fname), 'w') as f:
        f.write(param_log)

# -----------------

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
        c1 = int(math.ceil(0.5 * seqlen - 0.5 * corelen))
    spos = (s1, s1 + width)
    cpos = (c1, c1 + corelen)
    # process each core separately and make sure that the list is unique
    for core in set(corelist):
        seq_wcore = [(score, seq[spos[0]:spos[1]]) for score, seq in seqintensities if seq[cpos[0]:cpos[1]] == core]
        seq_core_df = pd.DataFrame(seq_wcore, columns = ['score', 'seq'])
        agg = seq_core_df.groupby(['seq'], as_index=False).median()
        core_dict[core] = agg[["score","seq"]].values.tolist()
    return core_dict

def benchmark_kfold(param_dict, rows, numfold, kmers=[1,2,3]):
    """
    run_kfold wrapper for benchmarking per fold
    """
    scc_list = [] # r squared
    mse_list = []

    kf = run_kfold(param_dict, rows, numfold, kmers=[1,2,3])

    for i in range(len(kf)):
        p_label, p_acc, p_val = kf[i]["svmpred"]
        acc, mse, scc = p_acc
        scc_list.append(scc)
        mse_list.append(mse)

    avg_scc = sum(scc_list)/numfold
    avg_mse = sum(mse_list)/numfold

    return {"params":param_dict, "avg_scc": avg_scc, "avg_mse":avg_mse}

def predict_kfold(param_dict, rows, numfold, kmers=[1,2,3]):
    """
    run_kfold wrapper for predictions per fold
    """
    dflist = []
    kf = run_kfold(param_dict, rows, numfold, kmers=[1,2,3])
    for i in range(len(kf)):
        predicted = kf[i]["svmpred"][0]
        seqs = [x[1] for x in kf[i]["test"]]
        measured = [x[0] for x in kf[i]["test"]]
        ldict = [{"sequence":seqs[i],"measured":measured[i],"predicted":predicted[i], "fold":i+1} for i in range(len(seqs))]
        dflist.extend(ldict)
    return dflist

def run_kfold(param_dict, rows, numfold, kmers=[1,2,3]):
    """
    Run k KFold

    Args:
        param_dict: dictionary mapping param string to its value
        rows: input rows
        numfold: k for cross validation
        kmers: list of kmers, default [1,2,3]
    Return:
        dictionary of model performance (SCC, MSE) if benchmark is True, else
        return predictions for each fold
    """
    kf = KFold(numfold, shuffle=True)
    splitted = kf.split(rows)
    param_str = "-s 3 -b 1 -q " # epsilon-SVR, prob estimate true, quiet mode
    param_str += " ".join(["-{} {}".format(k,v) for k,v in param_dict.items()])
    params = svmutil.svm_parameter(param_str)

    foldidx = 1
    fold_results = []
    for train_idx, test_idx in splitted:
        train_list = [rows[i] for i in train_idx]
        test_list = [rows[i] for i in test_idx]

        y_train, x_train = libsvm_generate_matrix(train_list, kmers)
        y_test, x_test = libsvm_generate_matrix(test_list, kmers)

        train_prob  = svmutil.svm_problem(y_train, x_train)

        model = svmutil.svm_train(train_prob, params)
        #svmutil.svm_save_model('model_name.model', m)
        # y is only needed when we need the model performance
        svmpred = svmutil.svm_predict(y_test, x_test, model, options="-q")
        fold_results.append({"test":test_list, "svmpred":svmpred})
    return fold_results

def test_param_comb(traindata, param_dict, numfold=10, kmers=[1,2,3]):
    param_log = ""
    run_params = {}
    for core in traindata:
        #print("Working for core: %s" % core)
        run_params[core] = benchmark_kfold(param_dict, rows=traindata[core], numfold=numfold, kmers=kmers)
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
    p_use = 0.999 if p == 1 else p # to avoid division by zero
    return  np.log(p_use/(1.0-p_use))

def genmodel_gridsearch(cores_centered, user_param, numworkers):
    # get all user params needed. We use this instead of direct input functions
    # to make integration with the SLURM version easier
    param_dict = user_param["grid"]
    numfold = user_param["numfold"]
    kmers = user_param["kmers"]

    # Make list of params for grid search
    params = list(param_dict.keys())
    combinations = itertools.product(*param_dict.values())
    combinations = [{params[i]:c[i] for i in range(len(c))} for c in combinations]

    # ----- RUNNING CROSS VALIDATION -----
    core_params = {}
    for core in cores_centered:
        print("Working for core: %s" % core)
        benchmark_kfold_partial = functools.partial(benchmark_kfold, rows=cores_centered[core], numfold=numfold, kmers=kmers)
        with cc.ProcessPoolExecutor(max_workers = numworkers) as executor:
            # TODO: update input to combinations to dictionary
            run_params = list(tqdm(executor.map(benchmark_kfold_partial, combinations), total=len(combinations)))
        core_params[core] = run_params
    write_result(core_params, cores_centered, user_param)


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
