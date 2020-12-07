import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import pickle

def process_arrfile(df):
    df[["Name","type"]] = df["Name"].str.rsplit('_',n=1,expand=True)
    df = df.melt(id_vars=["Name", "type", "Sequence"], var_name="ori", value_name="affinity")
    df = df[~df["ori"].str.contains("Median")]
    df[["ori","rep"]] = df["ori"].str.split("_",expand=True)
    return df

def assign_label(l1, l2):
    if l1 == "cooperative" and l2 == "cooperative":
        return "cooperative"
    elif l1 == "anticoop" and l2 == "anticoop":
        return "anticoop"
    elif l1 == "additive" and l2 == "additive":
        return "additive"
    else:
        return "fail_cutoff"

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    probe_path = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd"
    dflist = [pd.read_csv("%s/2.processed_gpr/20191004_258614510001_ETS1_550_5_1-4_alldata.txt"%probe_path,sep="\t"),
           pd.read_csv("%s/2.processed_gpr/20191030_258614510001_550_55_ch4_alldata.txt"%probe_path,sep="\t")]
    dflist = [df[df["Name"].str.contains("ets1",na=False)].sort_values(by=["Name"]) for df in dflist]
    arr.plot_chamber_corr(dflist[0], dflist[1], median=True,
                           path="all_log.png", log=True,
                           title="Chamber 1-4 correlation", xlab="Chamber1", ylab="Chamber4")

    neg = pd.read_csv("%s/2.processed_gpr/20191004_258614510001_ETS1_550_5_1-4_alldata.txt"%probe_path,sep="\t")
    neg = neg[neg["Name"].str.contains("NegativeCtrl",na=False)][["Name","Alexa488Adjusted"]]
    neg[["Name","rep","ori"]] = neg["Name"].str.rsplit('_',n=2,expand=True)
    cutoff = float(neg.groupby(["Name","ori"]).median().reset_index()[["Alexa488Adjusted"]].quantile(0.95))# confirm
    print("Negcutoff", cutoff)

    arr_df = pd.read_csv("%s/3.coop_array_files/coop_20191004_258614510001_ETS1_550_5_1-4_alldata.tsv"%probe_path,sep="\t")
    arr_df = arr_df[~arr_df['Name'].str.contains(r'dist|weak')]
    arr_df = process_arrfile(arr_df)
    seqdf = arr_df[arr_df["type"] == "wt"][["Name","Sequence"]]
    indiv,two = arr.make_replicas_permutation(arr_df)
    # pickle.dump( indiv, open( "indivsum.p", "wb" ) )
    # pickle.dump( two, open( "twosites.p", "wb" ) )
    # indiv = pickle.load( open( "indivsum.p", "rb" ) )
    # two = pickle.load( open( "twosites.p", "rb" ) )

    lbled = arr.label_replicas_permutation(indiv, two, arr_df, cutoff=cutoff, pcut=0.01)
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
    arr.plot_classified_labels(lbled_both, col1="indiv_median", col2="two_median", log=True,
                       xlab="log(m1-m3+m2-m3)", ylab="log(wt-m3)", path="labeled_log_both.png", title="Cooperative plot, both orientations")
    seqlbled = lbled_both.merge(seqdf, on="Name")[["Sequence","label"]]
    seqlbled.to_csv("seqlbled.csv", index=False)
