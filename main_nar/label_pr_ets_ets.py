import chip2probe.training_gen.arranalysis as arr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

def assign_label(l1, l2):
    # We restrict same label in both orientations
    if l1 == "cooperative" and l2 == "cooperative":
        return "cooperative"
    elif l1 == "anticooperative" and l2 == "anticooperative":
        return "anticooperative"
    elif l1 == "independent" and l2 == "independent":
        return "independent"
    else:
        return "fail_cutoff"

def plot_pval(df, path, ori):
    dfin = pd.DataFrame(df[(df["label_%s" % ori] != "below_cutoff") & (df["label_%s" % ori] != "anticooperative")])[["p_%s"%ori, "label"]]
    dfin['p_%s' % ori] = dfin['p_%s' % ori].apply(lambda x: "%.4f" % x)
    dfin_c = dfin.groupby('p_%s' % ori)['p_%s' % ori].count().reset_index(name="count")
    dfin_c = dfin_c.set_index('p_%s' % ori)
    plt.rcParams["figure.figsize"] = (5,10)
    dfin_c.plot.barh(rot=0)
    plt.savefig(path)

# label the probes for ets1 ets1
if __name__ == "__main__":
    neg_percentile = 0.95
    pcut = 0.0035 # pval cutoff for classification
    fdrcorr = True
    pd.set_option("display.max_columns",None)

    # for the original array
    df, negdf = pd.read_csv("input/probefiles/Ets1_Ets1_pr_clean.csv"), pd.read_csv("input/probefiles/Ets1_Ets1_neg_clean.csv")
    df = df[df["Name"] != "ets1_GM23338_seq16_0"] # err
    nmsqmap = df[(df["ori"] == "o1") & (df["type"] == "wt")][["Name","Sequence"]].drop_duplicates()

    dftype = df[["Name","intensity","type"]].groupby(["Name","type"]).mean().reset_index()

    dftype_med = dftype.groupby(["Name","type"]).median().reset_index().pivot(index='Name', columns='type')
    dftype_med.columns = dftype_med.columns.droplevel(0)
    dftype_med = dftype_med.reset_index().merge(nmsqmap, on="Name")

    import sys
    sys.exit(0)

    print("Number of distinct names %d" % df["Name"].nunique())

    cutoff = 170.2740377474845
    """
    # Get the negative control cutoff
    seqnmneg = negdf[["Name","Sequence"]].drop_duplicates(subset="Name")
    negdf = negdf[["Name","intensity"]].groupby("Name") \
            .median().reset_index() \
            .merge(seqnmneg, on="Name")[["Sequence","intensity"]]
    negdf.to_csv("negdf.csv",index=False)
    cutoff = float(negdf[["intensity"]].quantile(neg_percentile))
    print("Negative control cutoff", cutoff)
    """

    # Make permutation of m1+m2-m3 (i.e. indiv) and wt (i.e. two)
    # indiv,two = arr.make_replicas_permutation(df, affcol="intensity")
    indiv, two = pickle.load( open( "indiv.p", "rb" ) ), pickle.load( open( "two.p", "rb" ) )

    arr.permutdict2df(indiv).drop_duplicates().sort_values(["Name","ori"]).to_csv("ets1_ets1_indiv.csv", index=False)
    arr.permutdict2df(two).drop_duplicates().sort_values(["Name","ori"]).to_csv("ets1_ets1_two.csv", index=False)

    # Make permutation of wt and m1+m2-m3
    lbled = arr.label_replicas_permutation(indiv, two, df, cutoff=cutoff, pcut=pcut, fdrcor=False, affcol="intensity")
    # import pickle
    # lbled = pickle.load( open( "save.p", "rb" ) )
    df["intensity"] = np.log(df["intensity"])
    lbled_both = lbled["o1"].merge(lbled["o2"], on="Name", suffixes=("_o1", "_o2"))
    lbled_both["label"] = lbled_both.apply(lambda x: assign_label(x["label_o1"], x["label_o2"]), axis=1)
    plot_pval(lbled_both, "p_o2_dist.png", "o2")

    seqdf = df[(df["type"] == "wt") & (df["ori"] == "o1")][["Name","Sequence"]]
    seqlbled = lbled_both.merge(seqdf, on="Name")[["Name","Sequence","label","p_o1","label_o1","p_o2","label_o2"]].drop_duplicates()
    seqlbled.to_csv("ets_ets_seqlabeled.csv", index=False)
    dftype_med.merge(seqlbled[seqlbled["label"] != "fail_cutoff"][["Name","label"]]).sort_values("Name")[["Sequence","m1", "m2", "m3", "wt", "label"]].to_csv("m1m2m3wt.csv",index=False)
    print("Label count", seqlbled["label"].value_counts())

    # plot both result in one orientation only; only take independent and cooperative since we have little to no anticooperative
    filt = lbled_both[(lbled_both['label'] != "fail_cutoff") & (lbled_both['label'] != "anticooperative")]
    lbled_one_ori = lbled['o1'].merge(filt[["Name"]],on="Name")

    lbled_one_ori.to_csv("lbled_o1_selected.csv", index=False)
    lbled_one_ori = pd.read_csv("lbled_o1_selected.csv")
    # arr.plot_classified_labels(lbled_one_ori, col1="indiv_median", col2="two_median", plotnonsignif=False,
    #                    xlab="M1-M3+M2-M3", ylab="WT-M3", path="labeled_log_one_both.pdf", title="Cooperative vs independent binding of Ets1-Ets1",
    #                    labelnames=["cooperative","independent","anticooperative"], showlog=True)
