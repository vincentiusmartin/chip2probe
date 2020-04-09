
import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe")

import pandas as pd
import random

from trainingdata.training import Training
from trainingdata.dnashape import DNAShape

from main import *


def make_ets1_mutations(df, flen=2):
    #mutdict = {}
    allmut = []
    atgc = {'A','T','G','C'}
    for idx,row in df.iterrows():
        if row["distance"] < 6: # we don't have enough flank
            continue
        mutlist = []
        seq = row["sequence"]
        s1pos = row["site_wk_pos"] if row["site_wk_pos"] < row["site_str_pos"] else row["site_str_pos"]
        s2pos = row["site_wk_pos"] if row["site_wk_pos"] > row["site_str_pos"] else row["site_str_pos"]
        s1 = seq[s1pos-2:s1pos+2]
        s2 = seq[s2pos-2:s2pos+2]
        for s in [s1,s2]:
            # first do it for the core
            if s == "GGAA":
                mutlist.append(seq[:s1pos-2] + "GGAT" + seq[s1pos+2:])
            elif s == "GGAT":
                mutlist.append(seq[:s1pos-2] + "GGAA" + seq[s1pos+2:])
            elif s == "TTCC":
                mutlist.append(seq[:s1pos-2] + "ATCC" + seq[s1pos+2:])
            elif s == "ATCC":
                mutlist.append(seq[:s1pos-2] + "TTCC" + seq[s1pos+2:])
            else:
                print("undefined core: %s" % s)
        for sp in [s1pos,s2pos]:
            # then do it for the flanks
            pos = list(range(sp-2-flen,sp-2)) + list(range(sp+2,sp+2+flen))
            for p in pos:
                for nuc in atgc - set(seq[p]):
                    mutlist.append(seq[:p] + nuc + seq[p+1:])
        wt_tbl = [{"sequence":row["sequence"], "site_wk_pos":row["site_wk_pos"], "site_str_pos":row["site_str_pos"], "distance":row["distance"], "wtlabel":row['label']}]
        mut_tbl = [{"sequence":m, "site_wk_pos":row["site_wk_pos"], "site_str_pos":row["site_str_pos"], "distance":row["distance"], "wtlabel":row['label']} for m in mutlist]
        allmut += wt_tbl + random.sample(mut_tbl,2)
        if idx > 300:
            break
    return allmut

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
    return ypred,yprob

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


if __name__ == '__main__':
    trainingpath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/training_data/training_p01_adjusted.tsv"
    pd.set_option('display.max_columns', None)
    #df = pd.read_csv(trainingpath, sep="\t")
    #md = make_ets1_mutations(df)

    #pd.DataFrame(md).to_csv("mutlist.csv", index=False)
    df = pd.read_csv("mutlist.csv")

    test = Training(df,4)
    #test = Training(pd.DataFrame(md),4)
    ds = DNAShape("all_model_files/dnashape/0")

    xt = test.get_feature_all({
        "distance":{"type":"numerical"},
        "flankshape": {"ds":ds, "seqin":5, "smode":"positional"},
        "flankshape": {"ds":ds, "seqin":-3, "smode":"positional"},
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    })

    idx1 = df.index[df['sequence'] == "GTTTGATCCAGGAAATGGTGTCCTTCCTGTGGACCT"].tolist()[0]
    idx2 = df.index[df['sequence'] == "GTTTGATCCAGGAAACGGTGTCCTTCCTGTGGACCT"].tolist()[0]
    print(xt[idx1])
    print(xt[idx2])

    # p1label, p1prob = predict_all_ori(test,ds)
    # p2label, p2prob = predict_per_ori(test,ds)
    # seqlist = df['sequence'].values.tolist()
    # wtlabel = df['wtlabel'].map({'cooperative': 1, 'additive': 0}).values.tolist()
    # labels = ["sequence", "y_wt", "ypred_all", "yprob_all", "ypred_ori", "ypred_ori"]
    # res = list(zip(seqlist, wtlabel, p1label, p1prob, p2label, p2prob))
    # pd.DataFrame(res,columns=labels).to_csv("mutpred.csv")

    # p2 = model2.predict(x_test2)
    # pl2 = list(zip(seqlist,p2))
    # for i in range(len(pl2)):
    #     print(i,pl2[i])
