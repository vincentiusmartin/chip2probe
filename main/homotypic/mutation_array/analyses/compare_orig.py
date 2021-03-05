import os
import pandas as pd
from plotcor import read_chfile
import chip2probe.training_gen.arranalysis as arr
import matplotlib.pyplot as plt
import numpy as np

os.chdir("../../../..")

def read_chfile(path, key):
    df = pd.read_csv(path, sep="\t")[["Name", "ID","Sequence","Alexa488Adjusted"]]
    df = df[df["ID"].str.contains(key, na=False)]
    df["Sequence"] = df["Sequence"].str[:36]
    negdf = df[df["Name"].str.contains("NegativeCtrl", na=False)]
    return df, negdf

def assign_label(l1, l2):
    if l1 == "cooperative" and l2 == "cooperative":
        return "cooperative"
    elif l1 == "anticoop" and l2 == "anticoop":
        return "anticoop"
    elif l1 == "additive" and l2 == "additive":
        return "additive"
    else:
        return "fail_cutoff"

def label_one_two(indiv, two):
    print(indiv.sort_values(["id","ori"]))
    for ori in ["o1", "o2"]:
        curindiv = indiv[indiv["ori"] == ori]
        curtwo = two[two["ori"] == ori]
    indiv_dict = indiv[["id","ori","affinity"]].groupby(["id","ori"])['affinity'].apply(list).to_dict()
    print(indiv_dict)

if __name__ == "__main__":
    # pd.set_option("display.max_columns",None)
    # sl_orig = pd.read_csv("output/homotypic/training/seqlbled.csv").drop_duplicates()
    # sl_cust = pd.read_csv("output/homotypic/custom_sequences/30nM/01_wfdr/seqlbled.csv")[["Sequence","label"]].drop_duplicates()
    # sl = sl_orig.merge(sl_cust, on="Sequence", suffixes=("_orig", "_cust"))
    # print("Number of same predictions: %d/%d" % (sl[sl["label_orig"] == sl["label_cust"]].shape[0], sl.shape[0]))
    # sl.to_csv("compare_orig.csv", index=False)
    # seqdf = sl[["Sequence"]].reset_index()
    # seqdf['index'] = seqdf.index

    # indiv_orig, two_orig = pd.read_csv("output/homotypic/training/indiv.csv"),  pd.read_csv("output/homotypic/training/two.csv")
    # indiv_cust, two_cust = pd.read_csv("output/homotypic/custom_sequences/20nM/indiv.csv"),  pd.read_csv("output/homotypic/custom_sequences/20nM/two.csv")
    # the sequence is the wildtype, so we can just use it
    # origseqs = set(indiv_orig["Sequence"].unique().tolist())
    # custseqs = set(indiv_cust["Sequence"].unique().tolist())
    # inter = origseqs.intersection(custseqs)
    # inter_idx = list(range(len(inter)))
    # interseqs = list(zip(inter_idx, inter))
    # interdf = pd.DataFrame(interseqs, columns=["id","Sequence"])
    # indiv_orig, two_orig = indiv_orig.merge(interdf, on="Sequence").sort_values("id"), two_orig.merge(interdf, on="Sequence").sort_values("id")
    # indiv_cust, two_cust = indiv_cust.merge(interdf, on="Sequence").sort_values("id"), two_cust.merge(interdf, on="Sequence").sort_values("id")
    # arr.plot_ori_inconsistency(indiv_orig, two_orig, namecol="id", prefix_path="orig", log=True)
    # arr.plot_ori_inconsistency(indiv_cust, two_cust, namecol="id", prefix_path="cust", log=True)

    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata"
    df_orig, neg_orig = read_chfile("%s/191030_coop-PBM_Ets1_v1_2nd/2.processed_gpr/20191004_258614510001_ETS1_550_5_1-4_alldata.txt" % basepath, "Coop1Ets")
    # 201128_validation_array_ets1_v2_1/10nMEts1_alexa488_550_20_alldata.txt %s/210102_validation_array_ets1_v2_2/20nMEts1_alexa488_550_10_alldata.txt
    df_cust, neg_cust = read_chfile("%s/210102_validation_array_ets1_v2_2/30nMEts1_alexa488_550_10_alldata.txt" % basepath, "Coop2Ets")
    pd.set_option("display.max_columns",None)
    df_orig = df_orig[df_orig["Name"].str.contains("_m1_") | df_orig["Name"].str.contains("_m2_") | df_orig["Name"].str.contains("_wt_")]
    df_cust = df_cust[df_cust["Name"].str.contains("_m1_") | df_cust["Name"].str.contains("_m2_") | df_cust["Name"].str.contains("_wt_")]
    inter = set(df_orig["Sequence"].unique().tolist()).intersection(set(df_cust["Sequence"].unique().tolist()))
    inter_idx = list(range(len(inter)))
    interseqs = list(zip(inter_idx, inter))
    interdf = pd.DataFrame(interseqs, columns=["id","Sequence"])
    df_orig = df_orig.merge(interdf, on="Sequence") \
                    .sort_values("id")
    df_cust = df_cust.merge(interdf, on="Sequence").sort_values("id") \
                    .sort_values("id")

    ax = plt.axes()
    arr.plot_chamber_corr(df_orig, df_cust, median=True, ax=ax,
                        namecol="id", log=True, xlab="first array (10nM)", ylab="validation array (30nM)")

    df_orig_wt = df_orig[df_orig["Name"].str.contains("_wt_")]
    df_cust_wt = df_cust[df_cust["Name"].str.contains("_wt_")]
    dfx = df_orig_wt.groupby(["id","Sequence"])[["Alexa488Adjusted"]].median().reset_index()
    dfy = df_cust_wt.groupby(["id","Sequence"])[["Alexa488Adjusted"]].median().reset_index()
    dfcombined = dfx.merge(dfy, on=["id","Sequence"])[["id","Sequence","Alexa488Adjusted_x", "Alexa488Adjusted_y"]]
    dfcombined["Alexa488Adjusted_x"] = np.log(dfcombined["Alexa488Adjusted_x"])
    dfcombined["Alexa488Adjusted_y"] = np.log(dfcombined["Alexa488Adjusted_y"])

    selected = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/array_design_files/Coop2Ets_validation/custom_probes_selected.csv")
    wtseqs = selected[selected["comment"] == "wt"][["sequence","wtlabel"]].drop_duplicates().rename({"sequence":"Sequence"}, axis=1)

    for lab in [0,1]:
        wtlab = wtseqs[wtseqs['wtlabel'] == lab][['Sequence']]
        cmb = dfcombined.merge(wtlab)

        c, l = ("red", "wt_add") if lab == 0 else ("blue", "wt_coop")
        ax.scatter(cmb["Alexa488Adjusted_x"].values, cmb["Alexa488Adjusted_y"].values, s=3, color=c, label=l)
    ax.legend(loc="lower right")
    plt.savefig("cor_30nM_with_lbled_wt.png")
