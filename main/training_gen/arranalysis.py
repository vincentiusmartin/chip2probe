import os

import pandas as pd

import chip2probe.training_gen.arranalysis as arr
import matplotlib.pyplot as plt
import chip2probe.util.stats_r as st

# this is not a general function, so we put in main
def fix_naming(df):
    df[["Name","rep","ori"]] = df["Name"] \
        .str.rsplit("_", n = 2, expand = True)
    df = df.drop(columns=["ori"])

    df_ori = df[["Name","Sequence"]] \
        .drop_duplicates()
    df_ori["ori"] = df_ori \
        .groupby("Name",as_index=False) \
        .cumcount() + 1
    df_ori["ori"] = "o" + df_ori["ori"].astype(str)

    df = df.merge(df_ori, on=["Name","Sequence"])
    return df

if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/coop_hetero-PBM_Ets_EtsRunx_v1"
    dflist = [pd.read_csv("%s/Ets1_70.txt"%basepath,sep="\t"),
            pd.read_csv("%s/Ets1_Runx1_70.txt"%basepath,sep="\t")]

    # ------- plot negative control -------
    tf_str = "runx1"
    ori = "o1"
    neglist = []
    for df in dflist:
        df = df[df["Name"].str.contains("seq",na=False)] \
                .sort_values(by=["Name"])
        # fix the naming error
        negdf = pd.DataFrame(df[df["Name"].str.contains("negative_control")])
        negdf = fix_naming(negdf)
        negdf = pd.DataFrame(negdf[negdf["Name"].str.contains(tf_str)]) #  & (negdf["ori"] == ori)
        neglist.append(negdf)
    #arr.plot_chamber_corr(neglist[0], neglist[1], median=True, log=True , extrajoincols=["ori"])

    # get the negative control cutoff, we do it from first chamber
    cutoff = float(neglist[0][["Alexa488Adjusted"]].quantile(0.95))
    print(cutoff)

    # ------- labeling sequences, for now only take wt -------
    keyword = "all_clean_seqs"
    df_genomics = []
    for df in dflist:
        df_gen = df[df["Name"].str.contains(keyword,na=False)].sort_values(by=["Name"])
        df_gen[["Name","type","rep","ori"]] = df_gen["Name"].str.rsplit("_", n = 3, expand = True)
        df_gen = df_gen[df_gen["type"] == "wt"] #only take wt
        df_genomics.append(df_gen)

    for ori in ["o1","o2"]:
        df1 = df_genomics[0][df_genomics[0]["ori"] == ori]
        df2 = df_genomics[1][df_genomics[1]["ori"] == ori]
        break
