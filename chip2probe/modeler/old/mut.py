
import sys

import pandas as pd
import random
import math
import logging


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
            comment = "2_right_" + str(p - (spos[1] + 2) + 2)
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
    pass

if __name__ == '__main__':
    #-----------------generate list of mutated sequences------------------
    # trainingpath = "../../../input/modeler/training_data/training_p01_adjusted.tsv"
    # pd.set_option('display.max_columns', None)
    # df = pd.read_csv(trainingpath, sep="\t")

    modelpaths = ["../../../input/modeler/imads_model/Ets1_w12_GGAA.model",
                  "../../../input/modeler/imads_model/Ets1_w12_GGAT.model"]
    modelcores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(modelpath, modelcore, 12, [1, 2, 3])
                    for modelpath, modelcore in
                    zip(modelpaths, modelcores)]
    imads = iMADS(imads_models, imads_threshold=IMADS_THRESHOLD)
    # #print(imads.predict_sequence("GTTTGATCCAGGAAATGGTGTCCTTCCTGTGGACCT"))
    # md = make_ets1_mutations(df, imads, flanklen=4)
    # pd.DataFrame(md).to_csv("mutlist.csv", index=False)

    # --------------predict labels of new mutated sequences-----------------
    # we have 72833 rows, including wild types
    df = pd.read_csv("mutlist.csv")
    df = df[:100]
    t = Training(df,corelen=4)
    # #test = Training(pd.DataFrame(md),4)
    # ds = DNAShape("all_model_files/dnashape/0")
    #
    xt = t.get_feature_all({
        "distance":{"type":"numerical"},
        "sitepref": {"imadsmodel": imads, "modelwidth":12}, # use custom imads
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    })
    xt = pd.DataFrame(xt).values.tolist()
    #print(xt)
    
    model = pickle.load(open("dist_ori_12merimads.sav", "rb"))
    ypred = model.predict(xt)
    yprob = model.predict_proba(xt)
    yprob = [yprob[i][0] if ypred[i] == 0 else yprob[i][1] for i in range(len(ypred))]
    print(ypred)
    # p1label, p1prob = predict_all_ori(test,ds)
    # p2label, p2prob = predict_per_ori(test,ds)
    seqlist = df['sequence'].values.tolist()
    wtlabel = df['wtlabel'].map({'cooperative': 1, 'additive': 0}).values.tolist()
    labels = ["sequence", "y_wt", "ypred_all", "yprob_all"] #,  "ypred_ori", "ypred_ori"]
    res = list(zip(seqlist, wtlabel, ypred, yprob)) #p1label, p1prob, p2label, p2prob))
    pd.DataFrame(res,columns=labels).to_csv("mutpred.csv")

    # p2 = model2.predict(x_test2)
    # pl2 = list(zip(seqlist,p2))
    # for i in range(len(pl2)):
    #     print(i,pl2[i])
