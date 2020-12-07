import os
os.chdir("../..")

import pandas as pd

import chip2probe.training_gen.arranalysis as arr
import matplotlib.pyplot as plt
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
    seqdf.to_csv("seq_w1_runx_site.csv")
    return seqdf


#d1, d2 = pd.read_csv("%s/ets1_only_median_ch1.csv"%basepath), pd.read_csv("%s/ets1_only_median_ch2.csv"%basepath)
#arr.plot_chamber_corr(d1,d2,valcol="intensity",extrajoincols=["orientation"],title="Probes with Ets1 sites (using m1/m2)")

if __name__ == "__main__":
    params = {'axes.labelsize': 22,
          'axes.titlesize': 16,
          "xtick.labelsize" : 14, "ytick.labelsize" : 14 , "axes.labelsize" : 14}
    plt.rcParams.update(params)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/coop_hetero-PBM_Ets_EtsRunx_v1"
    dflist = [pd.read_csv("%s/Ets1_70.txt"%basepath,sep="\t"),
           pd.read_csv("%s/Ets1_Runx1_70.txt"%basepath,sep="\t")]
    # dflist = [pd.read_csv("%s/Runx1_80.txt"%basepath,sep="\t"),
    #         pd.read_csv("%s/Runx1_Ets1_80.txt"%basepath,sep="\t")]
    dflist = [df[df["Name"].str.contains("seq",na=False)].sort_values(by=["Name"]) for df in dflist]

    kompas_ets = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)
    kompas_runx = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Runx1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)

    # ------- plot negative control -------
    #ori = "o1"
    filteredlist = []

    # get sequence with the site we want, just need to run this once
    #seqdf = pd.DataFrame(dflist[0][dflist[0]["Name"].str.contains("m1|m2")])
    #get_seq_wsite(seqdf, kompas_runx).to_csv("seq_w1_runx_site.csv")
    tf_bound = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/seq_w1_ets_site.csv")
    # tf_bound = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/seq_w1_runx_site.csv")

    pd.set_option('display.max_columns', None)
    tf_str = "ets1"
    pattern = "m1|m2" #"m1|m2" # negative_control
    for i in range (0,len(dflist)):
        df = dflist[i]
        # filtdf = pd.DataFrame(df[df["Name"].str.contains(tf_str) & df["Name"].str.contains(pattern)]) #for negctrl
        # filtdf = fix_naming(filtdf) # to fix the naming error
        # if i == 1: # normalization
        #    filtdf["Alexa488Adjusted"] = (filtdf["Alexa488Adjusted"] - 108.83)/0.87 #(df["Alexa488Adjusted"] - 107.4) / 0.82
        filtdf = pd.DataFrame(df[df["Name"].str.contains(pattern)]) # for  non negctrl
        filtdf = fix_naming(filtdf).merge(tf_bound, on=["Sequence"]) # for  non negctrl
        filtdf = filtdf[filtdf["has_site"] == True].groupby(["Name","ori"]).median().reset_index() # for non negctrl
        filteredlist.append(filtdf)
    titleplt = "Probes with Runx sites (using m1/m2)" # Chamber3 (Runx) vs Chamber4 (Runx+Ets) in log val
    arr.plot_chamber_corr(filteredlist[0], filteredlist[1], median=True,
                           extrajoincols=["ori"], path="%s_log.png"%pattern, log=False,
                           title="Sequences with 1 Ets1 site in Ets1-Runx1", xlab="Ets1 + Ab_Ets1 (x10⁴)", ylab="Ets1 + Runx1 + Ab_Ets1 (x10⁴)",)

    # get the negative control cutoff, we do it from first chamber
    allmed_ch1 = filteredlist[0].groupby(["Name","ori"]).median().reset_index()
    cutoff = allmed_ch1[["Alexa488Adjusted"]].quantile(0.75)
    print(cutoff)
