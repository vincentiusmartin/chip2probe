import os
os.chdir("../../../..")

import chip2probe.training_gen.arranalysis as arr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    pd.set_option("display.max_columns", None)
    targetdf = pd.read_csv("output/homotypic/custom_sequences/20nM/05_nofdr_either_ori/wtmatch/matchwt_nopos_add.csv")
    targetdf = targetdf[targetdf["muttype"] == "distance"].drop_duplicates()
    print(targetdf)
    targetdf = targetdf.rename({"id":"Name"},axis=1)
    trgtids = targetdf[["Name"]].drop_duplicates().sort_values("Name")

    # df, neg = arr.read_chamber_file("../probedata/210102_validation_array_ets1_v2_2/20nMEts1_alexa488_550_10_alldata.txt", "Coop2Ets")
    # df = df.rename({"Alexa488Adjusted":"affinity"}, axis=1)
    # df = df.merge(trgtids,on="Name")
    # seqdf = df[(df["type"] == "wt") & (df["ori"] == "o1")][["Name","Sequence"]].drop_duplicates()
    # indiv,two = arr.make_replicas_permutation(df)
    # ind_df = arr.permutdict2df(indiv)
    # two_df = arr.permutdict2df(two) #.merge(seqdf, on="Name").drop_duplicates()
    # ind_mean = ind_df.groupby("Name").mean().reset_index()
    # two_mean = two_df.groupby("Name").mean().reset_index()
    # mean_delta = two_mean.merge(ind_mean, on="Name")
    # mean_delta["delta"] = mean_delta["affinity_x"] - mean_delta["affinity_y"]

    # idst = targetdf[["Name","distance"]]
    # idst = idst.merge(mean_delta, on="Name")
    # idst.plot(x='distance', y='delta', style='o')
    # plt.savefig("here.png")
