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

def read_chfile(path, key):
    df = pd.read_csv(path, sep="\t")[["Name", "ID","Sequence","Alexa488Adjusted"]]
    df = df[df["ID"].str.contains(key, na=False)]
    df = df[~df['Name'].str.contains(r'dist|weak')]
    df["Sequence"] = df["Sequence"].str[:36]
    negdf = df[df["Name"].str.contains("NegativeCtrl", na=False)]
    negdf[["Name","rep","ori"]] =  negdf["Name"].str.rsplit("_", n = 2, expand = True)
    df = df[~df["Name"].str.contains("NegativeCtrl", na=False)]
    df[["Name","type","rep","ori"]] = df["Name"].str.rsplit("_", n = 3, expand = True)
    return df, negdf

# val array: /Users/vincentiusmartin/Research/chip2gcPBM/probedata/201128_validation_array_ets1_v2_1
# orig_file: /Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd
if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/210102_validation_array_ets1_v2_2"
    df, neg = read_chfile("/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/2.processed_gpr/20191004_258614510001_ETS1_550_5_1-4_alldata.txt", "Coop1Ets")
    cutoff = float(neg.groupby(["Name","ori"]).median().reset_index()[["Alexa488Adjusted"]].quantile(0.95))# confirm
    print("Negcutoff", cutoff)
    df = df.rename({"Alexa488Adjusted":"affinity"}, axis=1)
    seqdf = df[(df["type"] == "wt") & (df["ori"] == "o1")][["Name","Sequence"]].drop_duplicates()
    idmap = df[(df["ori"] == "o1") & (df["rep"] == "r1")]
    idmap["id"] =  df.apply(lambda x: "%s_%s" % (x["Name"],x["type"]), axis=1)
    idmap = idmap[["id","Sequence"]]
    seqdf["Name"] = seqdf["Name"].astype('str')
    seqdf.sort_values("Name").to_csv("seqdf.csv",index=False)

    print(df.shape)
    indiv,two = arr.make_replicas_permutation(df)
    arr.permutdict2df(indiv).merge(seqdf, on="Name").drop_duplicates().to_csv("indiv.csv", index=False)
    arr.permutdict2df(two).merge(seqdf, on="Name").drop_duplicates().to_csv("two.csv", index=False)
    # pickle.dump( indiv, open( "indivsum.p", "wb" ) )
    # pickle.dump( two, open( "twosites.p", "wb" ) )
    # indiv = pickle.load( open( "indivsum.p", "rb" ) )
    # two = pickle.load( open( "twosites.p", "rb" ) )

    lbled = arr.label_replicas_permutation(indiv, two, df, cutoff=cutoff, fdrcor=True, pcut=0.01)
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
    print("Count both\n",lbled_both["label"].value_counts())
    arr.plot_classified_labels(lbled_both, col1="indiv_median", col2="two_median", log=True,
                       xlab="log(m1-m3+m2-m3)", ylab="log(wt-m3)", path="labeled_log_both.png", title="Cooperative plot, both orientations")
    lbled_both.to_csv("lbled_both.csv", index=False)
    lbled_both["Name"] = lbled_both["Name"].astype('str')
    seqlbled = lbled_both \
                .merge(seqdf, on="Name")[["Sequence","label"]] \
                .merge(idmap, on="Sequence")[["id","Sequence","label"]]
    seqlbled.drop_duplicates().to_csv("seqlbled.csv", index=False)
