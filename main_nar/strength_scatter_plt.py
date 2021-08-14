import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

def get_median_wt(df,train,tfname):
    df = df[df["type"] == "wt"]
    dfmed = df.groupby("Name").median()["intensity"].reset_index()
    train = dfmed.merge(train[["Name","distance","orientation","label"]], on="Name") \
                 .rename({"intensity":"%s_intensity"%tfname},axis=1)
    # train["runx_intensity"] = np.log(train["runx_intensity"])
    return train

pd.set_option("display.max_columns",None)
if __name__ == "__main__":
    df = pd.read_csv("output/Ets1Runx1/training/train_ets1_runx1.tsv", sep="\t")
    print(df)
    # runx_raw = pd.read_csv("input/probefiles/Runx1_only_pr_clean.csv")
    # df = get_median_wt(runx_raw,dftrain,"runx")

    # df = pd.read_csv("etc/withimads.csv")[["Name","runx1_score"]] \
    #         .merge(dftrain[["Name","distance","orientation","label"]], on="Name")

    labels = df["label"].drop_duplicates().tolist()
    targetcol = "runx1_score"
    minrow = 5
    lbltarget = "independent"

    coop_prop = df.groupby(["distance","orientation"]) \
        .apply(lambda g: pd.Series({"prop":g[g["label"] == lbltarget].shape[0]/g.shape[0]})) \
        .reset_index()

    coop_med = df[df["label"] == lbltarget].groupby(["distance","orientation"])[targetcol] \
        .agg(["count","median"]) \
        .reset_index()

    coop_df = coop_med[["distance","orientation","count","median"]].merge(coop_prop, on=["distance","orientation"])
    coop_df.to_csv("strength_%s.csv"%lbltarget,index=False,float_format='%.2f')

    coop_df = coop_df[coop_df["orientation"] == "+/+"]
    coop_df = coop_df[coop_df["count"] >= minrow]
    coop_df["lbl"] = coop_df.apply(lambda x: "d=%s,%s" % (x["distance"],x["orientation"]), axis=1)
    coop_df.plot.scatter(x='median', y='prop', c="black")

    # for i, row in coop_df.iterrows():
    #     plt.annotate(row["lbl"], (row["median"], row["prop"]), fontsize=10)

    x,y = coop_df["median"], coop_df["prop"]
    plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), color='black', linestyle='dashed')
    slope, intercept = np.polyfit(x, y, 1)
    print(slope,intercept)
    abline_values = [slope * i + intercept for i in x]

    # plt.plot((min(x),max(abline_values)), (max(x),min(abline_values)))
    plt.xlabel('Runx1 strength median for %s probes' % lbltarget)
    plt.ylabel('Fraction of %s probes' % lbltarget)

    slope, intercept, r, p, std_err = stats.linregress(x,y)
    print('R = %0.3f, p=%0.3f' % (r,p))
    sign = '+' if intercept > 0 else '-'
    print('best fit: y=%0.2fx %s %0.2f' % (slope,sign,abs(intercept)))
    # plt.text(0.8, 0.8, 'R = %0.2f\np = %0.2f' % (r,p), color='red', fontsize=12, transform = plt.gca().transAxes)
    # sign = '+' if intercept > 0 else '-'
    # plt.text(0.8, 0.9, 'best fit: y=%0.2fx %s %0.2f' % (slope,sign,abs(intercept)), color='red', fontsize=12, transform = plt.gca().transAxes)

    fig = plt.gcf()
    # fig.set_size_inches(8, 4)
    fig.savefig('strength_%s_nolabel.png' % lbltarget)
