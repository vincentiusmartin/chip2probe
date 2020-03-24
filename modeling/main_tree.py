
import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe")

import pandas as pd

from trainingdata.training import Training
from trainingdata.dnashape import DNAShape

from main import *

idxtest = 30

def make_ets1_mutations(df, flen=2):
    mutdict = {}
    atgc = {'A','T','G','C'}
    for idx,row in df.iterrows():
        if idx != idxtest:
            continue
        print(row)
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
                mutlist.append(seq[:s1pos-2] + "GGAT" + seq[s1pos+2:])
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
        mut_tbl = [{"sequence":m, "site_wk_pos":row["site_wk_pos"], "site_str_pos":row["site_str_pos"], "distance":row["distance"]} for m in mutlist]
        mutdict[idx] = mut_tbl
        break
    return mutdict


if __name__ == '__main__':
    trainingpath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/training_data/training_p01_adjusted.tsv"
    pd.set_option('display.max_columns', None)
    df = pd.read_csv(trainingpath, sep="\t")
    md = make_ets1_mutations(df)

    model1 = pickle.load(open("all_model_files/model1_all_all_dt.sav", "rb"))
    model2 = pickle.load(open("all_model_files/model2_all_all_rf.sav", "rb"))

    test = Training(pd.DataFrame(md[idxtest]),4)
    ds = DNAShape("all_model_files/dnashape/0")
    x_dist = test.get_feature_distance(type="numerical")
    x_ori = test.get_feature_orientation(["GGAA","GGAT"], one_hot = "True")

    x_shape_in = test.get_feature_flank_shapes(ds, seqin = 5, site_mode="strength")
    x_shape_out = test.get_feature_flank_shapes(ds, seqin = -4, site_mode="strength")
    x_flank_in = test.get_feature_flank_core(k=3, seqin = 5, site_mode="strength")
    x_flank_out = test.get_feature_flank_core(k=3, seqin = -4, site_mode="strength")

    seqlist = [x["sequence"] for x in md[idxtest]]

    xtr1 = []
    for x in [x_dist, x_ori, x_flank_in, x_flank_out]: #
        xtr1 = merge_listdict(xtr1, x)
    x_test1 = pd.DataFrame(xtr1).values.tolist()

    xtr2 = []
    for x in [x_dist, x_ori, x_shape_in, x_shape_out]: #
        xtr2 = merge_listdict(xtr2, x)
    x_test2 = pd.DataFrame(xtr2).values.tolist()


    p1 = model1.predict(x_test1)
    pl1 = list(zip(seqlist,p1))
    for i in range(len(pl1)):
        print(i,pl1[i])
    print("\n")
    p2 = model2.predict(x_test2)
    pl2 = list(zip(seqlist,p2))
    for i in range(len(pl2)):
        print(i,pl2[i])
