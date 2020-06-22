
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import logging
import pickle

from chip2probe.probe_generator.probefilter.sitespredict.imads import iMADS
from chip2probe.probe_generator.probefilter.sitespredict.imadsmodel import iMADSModel
from chip2probe.modeler.old.training.training import Training
#from main import *

"""
Selection criteria:
    1. change the core ( 2 mutants )
    2. different mutation, similar preference
    3.
"""

IMADS_THRESHOLD = 0.19
#IMADS_THRESHOLD = 0.2128

def get_similar_score(mutlist, idx, rel_tol=1e-2):
    sim_idx = []
    mt = [(key, mutlist[key][idx]) for key in mutlist]
    srt = sorted(mt, key=lambda k: k[1]['score'])
    justadded = False
    for i in range(1, len(srt)):
        if math.isclose(srt[i-1][1]["score"], srt[i][1]["score"], rel_tol=rel_tol):
            if not justadded:
                sim_idx.append(srt[i-1][0])
            sim_idx.append(srt[i][0])
            justadded = True
        else:
            justadded = False
    return sim_idx


def get_valid_mutants(muts, imads):
    """
    Get mutated sequences with 2 binding sites.

    return: a list of dictionaries of sequence, comment, predict sequence dicts
    """
    mlist = []
    # get imads predictions for the cores in each sequence
    mlist = {i: imads.predict_sequence(muts[i]["sequence"]) for i in range(len(muts))}
    # filter to those that have 2 sites
    ret = []
    for idx in mlist:
        if len(mlist[idx]) == 2:
            # get the sequence, the two predictions for the two cores, and comment
            d = {"sequence": muts[idx]["sequence"],
                 "core1": mlist[idx][0],
                 "core2": mlist[idx][1],
                 "comment": muts[idx]["comment"]}
            ret.append(d)
    
    return ret

# math.isclose(a, b, rel_tol=1e-5)
def similar_score(wt, muts, imads, mlist, rndcount=2):
    wtpred = imads.predict_sequence(wt)
    mlist = {i: imads.predict_sequence(muts[i]) for i in range(len(muts))}
    # filter to those that have 2 sites
    mlist = [mlist[k] for k in mlist if len(mlist[k]) == 2]
    s1_muts, s2_muts = {}, {}
    for i in range(len(mlist)):
        sc1, sc2 = mlist[i][0]["score"], mlist[i][1]["score"]
        wtsite = 0 if math.isclose(sc1, wtpred[0]["score"], rel_tol=1e-2) else 1
        mutsite = 0 if math.isclose(sc2, wtpred[1]["score"], rel_tol=1e-2) else 1
        if wtsite == mutsite:
            continue
        if mutsite == 0:
            s1_muts[i] = mlist[i]
        else:
            s2_muts[i] = mlist[i]

    # get s1 and s2 mutants that are close
    s1_idx = get_similar_score(s1_muts, 0, rel_tol=0.01)
    s1_muts = [muts[idx] for idx in s1_idx]
    s2_idx = get_similar_score(s2_muts, 1, rel_tol=0.01)
    s2_muts = [muts[idx] for idx in s2_idx]
    mlisct_cpy = list(muts)
    for elm in s1_muts + s2_muts:
        mlisct_cpy.remove(elm)
    picked = s1_muts + s2_muts
    if len(mlisct_cpy) > rndcount:
        picked.extend(random.sample(mlisct_cpy, rndcount))
    return picked


def make_ets1_mutations(df, imads, flanklen=4):
    """
    Generate ets1 sequences where the cores or flanking regions are mutated.

    df: Dataframe containing sequences of interest
    imads: imads object used for prediction
    flanklen: maximum flank to mutate from the beginning/end of a binding site

    return: a list of dictionaries of the wild type and mutated sequences
            with corresponding features
    """
    allmut = []

    for idx, row in df.iterrows():
        if idx % 10 == 0:
            print("Processed %d/%d rows" % (idx, len(df)))
        # if they share the same flank, continue
        if row["distance"] < 12:
            continue
        seq = row["sequence"]
        s1pos = row["site_wk_pos"] if row["site_wk_pos"] < row["site_str_pos"] else row["site_str_pos"]
        s2pos = row["site_wk_pos"] if row["site_wk_pos"] > row["site_str_pos"] else row["site_str_pos"]
        spos = [s1pos, s2pos]

        # add the wildtype
        allmut += [{"sequence": row["sequence"], "site_wk_pos":row["site_wk_pos"], \
                    "site_str_pos":row["site_str_pos"], "site_wk_score": row["site_wk_score"],
                    "site_str_score": row["site_str_score"], "distance":row["distance"], \
                    "comment":"wt"}]

        # get list of mutations where the cores or flanks are mutated
        mutlist = mutate_core(seq, spos) + mutate_flanks(seq, spos, flanklen)

        # get list of mutated sequences with 2 binding sites
        mutlist = get_valid_mutants(mutlist, imads)

        # add all the valid mutated sequences
        for mut in mutlist:
            mid1 = mut["core1"]["core_mid"]
            mid2 = mut["core2"]["core_mid"]
            score1 = mut["core1"]["score"]
            score2 = mut["core2"]["score"]
            if score1 > score2:
                site_str_pos = mid1
                site_str_score = score1
                site_wk_pos = mid2
                site_wk_score = score2
            else:
                site_str_pos = mid2
                site_str_score = score2
                site_wk_pos = mid1
                site_wk_score = score1

            allmut += [{"sequence": mut["sequence"], "site_wk_pos": site_wk_pos, \
                        "site_str_pos":site_str_pos, "site_wk_score": site_wk_score,
                        "site_str_score": site_str_score, "distance":row["distance"], \
                        "comment": mut["comment"]}]


    return allmut


def mutate_core(seq, spos):
    """Given a sequence and the positions of the sites, find sequences where cores are mutated."""
    mutlist = []
    s1 = seq[spos[0] - 2:spos[0] + 2]
    s2 = seq[spos[1] - 2:spos[1] + 2]
    s = [s1, s2]

    comment = ""

    for i in [0, 1]:
        # add comment for the mutation
        if i == 0:
            comment = "core_1"
        else:
            comment = "core_2"

        if s[i] == "GGAA":
            newseq = seq[:spos[i]-2] + "GGAT" + seq[spos[i]+2:]
        elif s[i] == "GGAT":
            newseq = seq[:spos[i]-2] + "GGAA" + seq[spos[i]+2:]
        elif s[i] == "TTCC":
            newseq = seq[:spos[i]-2] + "ATCC" + seq[spos[i]+2:]
        elif s[i] == "ATCC":
            newseq = seq[:spos[i]-2] + "TTCC" + seq[spos[i]+2:]
        else:
            logging.warning("undefined core: %s" % s)
        # predict the new muatated sequence
        prednew = imads.predict_sequence(newseq)
        # get the sequence if both sites are sill present
        if len(prednew) == 2:
            mutlist.append({"sequence":newseq, 
                            "comment": comment})

    return mutlist


def mutate_flanks(seq, spos, flanklen):
    """Generated sequences with mutated flanking regions."""
    atgc = {'A', 'T', 'G', 'C'}
    mutflank = []
    # get the positions to mutate around the first binding site
    flank = min(spos[0] - 2, flanklen)
    pos = list(range(spos[0] - 2 - flank, spos[0] - 2)) + \
        list(range(spos[0] + 2, spos[0] + 2 + flanklen))
    # get the positions to mutate around the second binding site
    flank = min(len(seq) - spos[1] - 2, flanklen)
    pos += list(range(spos[1] - 2 - flanklen, spos[1] - 2)) + \
        list(range(spos[1] + 2, spos[1] + 2 + flank))
    # generate sequences with flanking regions mutated
    for p in pos:
        if spos[0]- 2 -flanklen <= p < spos[0]-2:
            comment = "1_left_" + str(spos[0] - 2 - p)
        elif spos[0] + 2 <= p < spos[0] + 2 + flanklen:
            comment = "1_right_" + str(p - (spos[0] + 2) + 1)
        elif spos[1] - 2 - flanklen <= p < spos[1] - 2:
            comment = "2_left_" + str(spos[1] - 2 - p)
        else:
            comment = "2_right_" + str(p - (spos[1] + 2) + 1)
        for nuc in atgc - set(seq[p]):
            mutflank.append({"sequence":seq[:p] + nuc + seq[p + 1:], 
                            "comment": comment})

    return mutflank

def predict_per_ori(test,ds):
    top_ht = ['dist_numeric','mgw_outer_wk_pos_0','mgw_outer_wk_pos_1','prot_outer_wk_pos_0', 'roll_outer_str_pos_2', 'prot_outer_wk_pos_2', 'helt_outer_str_pos_0', 'mgw_outer_wk_pos_2', 'helt_outer_str_pos_1' , 'helt_outer_wk_pos_2']
    top_hh =  [ 'dist_numeric', 'roll_outer_wk_pos_2', 'helt_outer_wk_pos_2', 'roll_outer_wk_pos_1', 'helt_outer_wk_pos_0', 'roll_outer_wk_pos_0', 'roll_outer_str_pos_0', 'helt_outer_wk_pos_1', 'prot_outer_str_pos_0', 'mgw_outer_str_pos_1']
    top_tt = ['dist_numeric','prot_outer_str_pos_0', 'mgw_outer_str_pos_0', 'roll_outer_wk_pos_2', 'mgw_outer_wk_pos_0', 'mgw_outer_wk_pos_2', 'prot_outer_wk_pos_2', 'helt_outer_str_pos_2', 'mgw_outer_wk_pos_1', 'mgw_outer_str_pos_1']
    model_ht = pickle.load(open("all_model_files/dist_ori_flank_ht.sav", "rb"))
    model_hh = pickle.load(open("all_model_files/dist_ori_flank_hh.sav", "rb"))
    model_tt = pickle.load(open("all_model_files/dist_ori_flank_tt.sav", "rb"))
    xt = test.get_feature_all({
        "distance":{"type":"numerical"},
        "flankshape": {"ds":ds, "seqin":5, "smode":"positional"},
        "flankshape": {"ds":ds, "seqin":-3, "smode":"positional"},
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":False}
    })
    ypred = []
    yprob = []
    for x_in in xt:
        if x_in['ori'] == "HH":
            xlist = [[x_in[k] for k in top_hh]]
            yp = model_hh.predict(xlist)
            pr = model_hh.predict_proba(xlist)[0]
        elif x_in['ori'] == "HT/TH":
            xlist = [[x_in[k] for k in top_ht]]
            yp = model_ht.predict(xlist)
            pr = model_ht.predict_proba(xlist)[0]
        elif x_in['ori'] == "TT":
            xlist = [[x_in[k] for k in top_tt]]
            yp = model_tt.predict(xlist)
            pr = model_tt.predict_proba(xlist)[0]
        ypred.append(yp[0])
        pr = pr[0] if yp == 0 else pr[1]
        yprob.append(pr)
    return ypred, yprob

def predict_all_ori(test,ds):
    model1 = pickle.load(open("all_model_files/dist_ori_flank_all.sav", "rb"))
    xt = test.get_feature_all({
        "distance":{"type":"numerical"},
        "flankshape": {"ds":ds, "seqin":5, "smode":"positional"},
        "flankshape": {"ds":ds, "seqin":-3, "smode":"positional"},
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    })
    xt_top10 = pd.DataFrame(xt)[['dist_numeric', 'prot_outer_wk_pos_1', 'ori_HH', 'prot_outer_str_pos_1', 'roll_outer_wk_pos_0', 'prot_outer_wk_pos_0', 'mgw_outer_wk_pos_2', 'helt_outer_wk_pos_2', 'prot_outer_str_pos_0', 'helt_outer_str_pos_1']]
    xt_top10 = xt_top10.values.tolist()
    ypred = model1.predict(xt_top10)
    yprob = model1.predict_proba(xt_top10)
    yprob = [yprob[i][0] if ypred[i] == 0 else yprob[i][1] for i in range(len(ypred))]
    return ypred,yprob

def populate_wt_pred(df):
    """Populate the dataframe with column indicating the wt pred."""
    df["wt_pred"] = 0
    for idx, row in df.iterrows():
        if row["comment"] == "wt":
            wt_pred = row["y_pred"]
        df.at[idx, "wt_pred"] = wt_pred

    return df

def filter_changed_pred_label(df):
    """Get rows in the df where wt pred and y pred are different."""
    df = df.loc[(df['wt_pred'] != df['y_pred']) | (df['comment']=='wt')]
    return df

def filter_high_pred_prob(df, cutoff=0.7):
    """Get rows in the df where y proba is at least cutoff."""
    df = df.loc[df['y_proba'] >= cutoff]
    return df

def filter_same_pred_label(df):
    """Get rows in the df where wt pred and y pred are the same."""
    df = df.loc[(df['wt_pred'] == df['y_pred'])]
    return df

def filter_similar_score(df, imads):
    """Get rows in the df where the imads scores for both sites are the same after mutation."""
    rows = []
    for idx, row in df.iterrows():
        # get the scores of the wild type sequence
        if row["comment"] == "wt":
            wt_seq = row["sequence"]
            wtpred = imads.predict_sequence(wt_seq)
        else:
            mut_seq = row["sequence"]
            seqpred = imads.predict_sequence(mut_seq)
            score1, score2 = seqpred[0]["score"], seqpred[1]["score"]
            if not (math.isclose(score1, wtpred[0]["score"], rel_tol=1e-2) and \
                math.isclose(score2, wtpred[1]["score"], rel_tol=1e-2)):
                rows.append(idx)

    # drop the rows where scores of corresponding sites are not similar
    df = df.drop(rows)
    return df
                
def clean_df(df):
    """Remove wild types with no eligible mutated sequences after filtering."""
    rows = []
    # number of mutated sequences for a particular wild type
    mut_count = 0
    wt_idx = 0
    for idx, row in df.iterrows():
        # get the scores of the wild type sequence
        if row["comment"] == "wt":
            # remove wild types with no eligible mutated sequence
            if idx > 0 and mut_count == 0:
                rows.append(wt_idx)
            wt_idx = idx
            mut_count = 0
        else:
            mut_count += 1

    df = df.drop(rows)
    return df

if __name__ == '__main__':
    
    trainingpath = "../../../input/modeler/training_data/training_p01_adjusted.tsv"
    modelpaths = ["../../../input/modeler/imads_model/Ets1_w12_GGAA.model",
                  "../../../input/modeler/imads_model/Ets1_w12_GGAT.model"]
    modelcores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(modelpath, modelcore, 12, [1, 2, 3])
                    for modelpath, modelcore in
                    zip(modelpaths, modelcores)]
    imads = iMADS(imads_models, imads_threshold=IMADS_THRESHOLD)

    #-----------------generate list of mutated sequences------------------
    pd.set_option('display.max_columns', None)
    df = pd.read_csv(trainingpath, sep="\t")

    md = make_ets1_mutations(df, imads, flanklen=4)
    pd.DataFrame(md).to_csv("mutlist.csv", index=False)

    # --------------predict labels of new mutated sequences-----------------
    #we have 72833 rows, including wild types
    df = pd.read_csv("mutlist.csv")
    
    t = Training(df,corelen=4)

    xt = t.get_feature_all({
        "distance":{"type":"numerical"},
        "sitepref": {"imadsmodel": imads, "modelwidth":12}, # use custom imads
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    })
    xt = pd.DataFrame(xt).values.tolist()
    
    model = pickle.load(open("dist_ori_12merimads.sav", "rb"))
    
    ypred = model.predict(xt)
    yprob = model.predict_proba(xt)
    yprob = [yprob[i][0] if ypred[i] == 0 else yprob[i][1] for i in range(len(ypred))]
    df["y_pred"] = ypred
    df["y_proba"] = yprob
    df = populate_wt_pred(df)
    df.to_csv("mutpred.csv", index=None)

    # --------Filter to get mutated sequences with changed label and high pred prob--------
    df = pd.read_csv("mutpred.csv")
    plt.hist(df["y_proba"], density=False, bins=10)
    plt.ylabel('Count')
    plt.xlabel('Probability')
    plt.title("Predicted Probability Distribution")
    plt.savefig("pred_prob_distribution.png")
    df = filter_changed_pred_label(df)
    df = filter_high_pred_prob(df, cutoff=np.percentile(df["y_proba"], 75))
    df.to_csv("mut_changed_label_high_pred_prob.csv", index=None)

    # ---------Filter to get mutated sequences where predicted labels don't change-------
    df = pd.read_csv("mutpred.csv")
    df = filter_same_pred_label(df)
    df.to_csv("mut_same_label.csv", index=None)

    #---------Filter to get sequences that have different predictions -------------
    #---------but have similar imads score for each corresponding site-------------
    df = pd.read_csv("mutpred.csv")
    df = filter_changed_pred_label(df)
    df = filter_similar_score(df, imads)
    df = clean_df(df)
    df.to_csv("mut_changed_label_similar_score.csv", index=None)


