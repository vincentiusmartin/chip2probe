import chip2probe.training_gen.arranalysis as arr
import pandas as pd
import numpy as np

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

# label the probes for ets1 ets1
if __name__ == "__main__":
    neg_percentile = 0.95
    pcut = 0.01 # pval cutoff for classification
    fdrcorr = True

    pd.set_option("display.max_columns",None)
    df, negdf = pd.read_csv("input/probefiles/Ets1_Ets1_pr_clean.csv"), pd.read_csv("input/probefiles/Ets1_Ets1_neg_clean.csv")
    df = df[df["Name"] != "ets1_GM23338_seq16_0"] # err
    print("Number of distinct names %d" % df["Name"].nunique())

    # Get the negative control cutoff
    cutoff = float(negdf.groupby(["Name","ori"]).median().reset_index()[["intensity"]].quantile(neg_percentile))
    print("Negative control cutoff", cutoff)

    # Make permutation of m1+m2-m3 (i.e. indiv) and wt (i.e. two)
    indiv,two = arr.make_replicas_permutation(df, affcol="intensity")

    # Make permutation of wt and m1+m2-m3
    lbled = arr.label_replicas_permutation(indiv, two, df, cutoff=cutoff, pcut=0.01, fdrcor=True, affcol="intensity")
    df["intensity"] = np.log(df["intensity"])
    lbled_both = lbled["o1"][["Name", "label", "indiv_median", "two_median"]].merge(lbled["o2"][["Name", "label", "indiv_median", 'two_median']], on="Name", suffixes=("_o1", "_o2"))
    lbled_both["label"] = lbled_both.apply(lambda x: assign_label(x["label_o1"], x["label_o2"]), axis=1)

    seqdf = df[(df["type"] == "wt") & (df["ori"] == "o1")][["Name","Sequence"]]
    seqlbled = lbled_both.merge(seqdf, on="Name")[["Name","Sequence","label"]].drop_duplicates()
    seqlbled.to_csv("ets_ets_seqlabeled.csv", index=False)
    print("Label count", seqlbled["label"].value_counts())

    # plot both result in one orientation only; only take independent and cooperative since we have little to no anticooperative
    filt = lbled_both[(lbled_both['label'] != "fail_cutoff") & (lbled_both['label'] != "anticooperative")]
    lbled_one_ori = lbled['o1'].merge(filt[["Name"]],on="Name")
    arr.plot_classified_labels(lbled_one_ori, col1="indiv_median", col2="two_median", log=True, plotnonsignif=False,
                       xlab="M1-M3+M2-M3", ylab="WT-M3", path="labeled_log_one_both.png", title="Cooperative vs independent binding of Ets1-Ets1",
                       labelnames=["cooperative","independent","anticooperative"])
    lbled_one_ori.to_csv("lbled_o1_selected.csv", index=False)
