import os
os.chdir("../..")

import pandas as pd

if __name__=="__main__":

    df_arr = pd.read_csv("Coop2Ets.tsv",header=None,sep="\t")[[1,2]]
    df_arr.columns = ["sequence_arr","id"]

    df_neg = df_arr.loc[df_arr["id"].str.contains("NegativeCtrl") &  df_arr["id"].str.contains("_o1_r1")]
    print("Num negctrl",len(df_neg["sequence_arr"].unique()))

    df_arr = df_arr.loc[df_arr["id"].str.contains("wt_o1_r1")]
    df_arr["id"] = [x[0] for x in df_arr["id"].str.rsplit('_',2).tolist()]

    print("Num wt sequence in the array",len(df_arr["sequence_arr"].unique()))
    orig_df = pd.read_csv("arrayseqs.csv")
    orig_df = orig_df.loc[orig_df["id"].str.contains("_wt")]
    orig = df_arr.merge(orig_df, on=["id"])[["id","sequence"]]
    print("Num wt sequence after first join",len(orig["sequence"].unique()))

    selected = pd.read_csv("custom_probes_selected.csv")
    print("Num wt sequence selected from the original pool", selected.shape[0])
    orig_back = orig.merge(selected, on=["sequence"])
    print("Num wt sequence after joining array with selected", len(orig_back["sequence"].unique()))
    print(orig_back)

    custom_count = orig_back.groupby(["select","muttype"])["seqid"] \
        .count()
    print("----------Number of sequences per group----------")
    print(custom_count)

    wt_count = orig_back.loc[orig_back["comment"] == "wt"]\
        .groupby(["select","muttype"])["seqid"] \
        .count()
    print("----------Number of wt per group----------")
    print(wt_count)
