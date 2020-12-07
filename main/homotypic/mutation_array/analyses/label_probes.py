from analyses import read_chfile
import chip2probe.training_gen.arranalysis as arr
import pickle

import pandas as pd

def assign_label(l1, l2):
    if l1 == "cooperative" and l2 == "cooperative":
        return "cooperative"
    elif l1 == "anticoop" and l2 == "anticoop":
        return "anticoop"
    elif l1 == "additive" and l2 == "additive":
        return "additive"
    else:
        return "fail_cutoff"

# val array: /Users/vincentiusmartin/Research/chip2gcPBM/probedata/201128_validation_array_ets1_v2_1
# orig_file: /Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd
if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/201128_validation_array_ets1_v2_1"
    cutoff = 120.99
    df, neg = read_chfile("%s/10nMEts1_alexa488_550_20_alldata.txt" % basepath)
    df = df.rename({"Alexa488Adjusted":"affinity"}, axis=1)
    seqdf = df[df["type"] == "wt"][["Name","Sequence"]].drop_duplicates()
    seqdf["Name"] = seqdf["Name"].astype('str')
    seqdf.to_csv("seqdf.csv",index=False)

    indiv,two = arr.make_replicas_permutation(df)
    # pickle.dump( indiv, open( "indivsum.p", "wb" ) )
    # pickle.dump( two, open( "twosites.p", "wb" ) )
    # indiv = pickle.load( open( "indivsum.p", "rb" ) )
    # two = pickle.load( open( "twosites.p", "rb" ) )

    lbled = arr.label_replicas_permutation(indiv, two, df, cutoff=cutoff, pcut=0.01)
    for ori in lbled:
        lbled[ori].to_csv("lbled_%s.csv"%ori,index=False)
        arr.plot_classified_labels(lbled[ori], col1="indiv_median", col2="two_median", log=True,
                        xlab="log(m1-m3+m2-m3)", ylab="log(wt-m3)", path="labeled_log_%s.png" % ori, title="Cooperative plot, %s" % ori)
        print("Count %s" % ori,lbled[ori]["label"].value_counts())

    lbled_both = lbled["o1"][["Name", "label", "indiv_median", "two_median"]].merge(lbled["o2"][["Name", "label", "indiv_median", 'two_median']], on="Name", suffixes=("_o1", "_o2"))
    lbled_both.to_csv("lbled_both.csv", index=False)
    lbled_both = pd.read_csv("lbled_both.csv")
    lbled_both["indiv_median"] = (lbled_both["indiv_median_o1"] + lbled_both["indiv_median_o2"]) / 2
    lbled_both["two_median"] = (lbled_both["two_median_o1"] + lbled_both["two_median_o2"]) / 2
    lbled_both["label"] = lbled_both.apply(lambda x: assign_label(x["label_o1"], x["label_o2"]), axis=1)
    print("Count",lbled_both["label"].value_counts())
    arr.plot_classified_labels(lbled_both, col1="indiv_median", col2="two_median", log=True,
                       xlab="log(m1-m3+m2-m3)", ylab="log(wt-m3)", path="labeled_log_both.png", title="Cooperative plot, both orientations")
    lbled_both.to_csv("lbled_both.csv", index=False)
    lbled_both["Name"] = lbled_both["Name"].astype('str')
    seqlbled = lbled_both.merge(seqdf, on="Name")[["Sequence","label"]]
    seqlbled.to_csv("seqlbled.csv", index=False)
