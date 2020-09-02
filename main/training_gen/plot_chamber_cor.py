import os
os.chdir("../..")

import pandas as pd

import chip2probe.training_gen.arranalysis as arr
import matplotlib.pyplot as plt
import chip2probe.util.stats_r as st
import statsmodels.stats.multitest as sm
from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

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

def get_seq_wsite(df, predictor):
    seqdf = df[["Sequence"]].drop_duplicates()
    seqdf["has_site"] = seqdf.apply(lambda x: len(predictor.predict_sequence(x["Sequence"])) == 1, axis=1)
    seqdf.to_csv("seq_w_site.csv")
    return seqdf


if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/coop_hetero-PBM_Ets_EtsRunx_v1"
    faricadir = ""
    dflist = [pd.read_csv("%s/Ets1_70.txt"%basepath,sep="\t"),
            pd.read_csv("%s/Ets1_Runx1_70.txt"%basepath,sep="\t")]

    kompas = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)

    # ------- plot negative control -------
    tf_str = "runx1"
    #ori = "o1"
    filteredlist = []
    pattern = "negative_control" #"m1|m2" # negative_control,

    # get sequence with the site we want, just need to run this once
    # get_seq_wsite(filtdf, kompas)
    tf_bound = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/seq_w1_ets_site.csv")

    for i in range (0,len(dflist)):
        df = dflist[i]
        df = df[df["Name"].str.contains("seq",na=False)] \
                .sort_values(by=["Name"])
        if i == 1:
            df["Alexa488Adjusted"] = (df["Alexa488Adjusted"] - 107.4) / 0.82
        # fix the naming error
        filtdf = pd.DataFrame(df[df["Name"].str.contains(tf_str) & df["Name"].str.contains(pattern)]) # for negctrl
        #filtdf = pd.DataFrame(df[df["Name"].str.contains(pattern)])
        filtdf = fix_naming(filtdf) #.merge(tf_bound, on=["Sequence"])
        #filtdf = filtdf[filtdf["has_site"] == True]
        filteredlist.append(filtdf)
    arr.plot_chamber_corr(filteredlist[0], filteredlist[1], median=True, extrajoincols=["ori"], path="negctrl_normalized.png", log=False, title="Ets1 negative control (normalized)")

    # get the negative control cutoff, we do it from first chamber
    cutoff = float(filteredlist[0][["Alexa488Adjusted"]].quantile(0.95))
    print(cutoff)
