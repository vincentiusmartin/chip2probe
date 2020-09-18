import os
os.chdir("../..")
import pandas as pd
import itertools
import matplotlib.pyplot as plt

from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel


if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/labeled/ch1_vs_ch2"
    df1 = pd.read_csv("%s/probes_labeled_er.csv"%basepath)[["Name","p_coop","p_anti","ets_pos","runx_pos","ori","label"]] \
        .drop_duplicates(ignore_index=True)
    df2 = pd.read_csv("%s/probes_labeled_re.csv"%basepath)[["Name","p_coop","p_anti","ets_pos","runx_pos","ori","label"]] \
        .drop_duplicates(ignore_index=True)
    df = df1.merge(df2, on=["Name"], suffixes=("_er", "_re"))

    # get probes that are additive in one orientation
    df_add1 = df[(df["label_er"] == "additive") & (df["label_re"] == "cooperative")][["p_coop_er"]].rename(columns={"p_coop_er":"p_coop"})
    df_add1["type"] = "additive_er"
    df_add2 = df[(df["label_er"] == "cooperative") & (df["label_re"] == "additive")][["p_coop_re"]].rename(columns={"p_coop_re":"p_coop"})
    df_add2["type"] = "additive_re"
    curdf = pd.concat([df_add1, df_add2])
    curdf["p_coop"] = curdf["p_coop"].astype(str)
    print(curdf)

    curdf.groupby('type')["p_coop"].value_counts().unstack(0).plot.barh()
    plt.title("p-val for additive in 1 orienttion only")
    plt.savefig("additive.png")

    # speed = [0.1, 17.5, 40, 48, 52, 69, 88]
    # index = ['snail', 'pig', 'elephant',
    #      'rabbit', 'giraffe', 'coyote', 'horse']
    # df = pd.DataFrame({'speed': speed}, index=index)
    # df.plot.bar(rot=0)
    # plt.savefig("additive.png")
