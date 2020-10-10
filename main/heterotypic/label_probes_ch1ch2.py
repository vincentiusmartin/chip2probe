import pandas as pd
import chip2probe.util.stats_r as st

from chip2probe.sitespredict.kompas import Kompas
import chip2probe.training_gen.arranalysis as arr
import matplotlib.pyplot as plt


def get_single_site(seq, predictor):
    s = predictor.predict_sequence(seq)
    if len(s) == 1:
        return s[0]['core_mid']
    else:
        return None

def get_sites_pos(df, pred1, pred2, tf1, tf2, o1name, o2name, seqcol="Sequence"):
    """
    Get site position for each sequence

    Args:
        df: input data frame

    """

    seqdf = df[[seqcol]].drop_duplicates()
    seqdf["%s_pos"%tf1] = seqdf.apply(lambda x: get_single_site(x["Sequence"], pred1),axis=1)
    seqdf["%s_pos"%tf2] = seqdf.apply(lambda x: get_single_site(x["Sequence"], pred2),axis=1)
    seqdf["ori"] = seqdf.apply(lambda x: o1name if x["%s_pos"%tf1] < x["%s_pos"%tf2] else o2name,axis=1)
    seqdf = seqdf[~seqdf['%s_pos'%tf1].isna() & ~seqdf['%s_pos'%tf2].isna()]
    return seqdf

def read_probe_data(df):
    df_gen = df[df["Name"].str.contains(keyword,na=False)].sort_values(by=["Name"])
    df_gen[["Name","type","rep","ori"]] = df_gen["Name"].str.rsplit("_", n = 3, expand = True)
    df_gen["Sequence"] = df_gen["Sequence"].str[0:36]
    df_gen = df_gen[df_gen["type"] == "wt"] #only take wt
    del df_gen["ori"]
    df_pos = get_sites_pos(df_gen, kompas_ets, kompas_runx, 'ets', 'runx', 'er', 're')
    df_probe = df_pos.merge(df_gen, on="Sequence")[["Name","Sequence","Alexa488Adjusted","ets_pos","runx_pos","ori","rep"]]
    return df_probe

if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/coop_hetero-PBM_Ets_EtsRunx_v1"
    # dflist = [pd.read_csv("%s/Runx1_80.txt"%basepath,sep="\t"),
    #        pd.read_csv("%s/Runx1_Ets1_80.txt"%basepath,sep="\t")]
    dflist = [pd.read_csv("%s/Ets1_70.txt"%basepath,sep="\t"),
            pd.read_csv("%s/Ets1_Runx1_70.txt"%basepath,sep="\t")]
    cutoff = 568.11 # from plot_chamber_corr.py
    ch_x = "Chamber1"
    ch_y = "Chamber2"
    p_default = 0.061
    p_ambiguous = 0.11

    kompas_ets = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)
    kompas_runx = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Runx1_kmer_alignment.txt",
                    core_start = 12, core_end = 17, core_center = 14)

    # ------- labeling sequences, for now only take wt -------
    keyword = "all_clean_seqs"
    df_genomics = []
    # print("Reading probe data")
    # for i in range(len(dflist)):
    #     df = dflist[i]
    #     df_gen = read_probe_data(df)
    #     if i == 1:
    #       df_gen["Alexa488Adjusted"] = (df_gen["Alexa488Adjusted"] - 108.83) / 0.87 # normalize
    #     df_gen.to_csv("df_wt_ch%d.csv" % (i+1), index=False)
    #     df_gen["Sequence"] = df_gen["Sequence"].str[0:36]
    #     df_genomics.append(df_gen)
    df_genomics = [pd.read_csv("df_wt_ch1.csv").convert_dtypes(), pd.read_csv("df_wt_ch2.csv").convert_dtypes()]

    print("Filtering")
    # we take sequence if either is above cutoff
    olist = []
    for ori in ["er","re"]:
        print("checking orientation %s" % ori)
        df1 = df_genomics[0][df_genomics[0]["ori"] == ori][["Name","rep","Alexa488Adjusted",]]
        df2 = df_genomics[1][df_genomics[1]["ori"] == ori][["Name","rep","Alexa488Adjusted"]]
        df_comb = df1.merge(df2, on=["Name","rep"])
        print("Number of distinct names before filtering %d" % df_comb["Name"].nunique())
        median_df = df_comb.groupby(["Name"]).median().reset_index()
        above_cut = median_df[(median_df["Alexa488Adjusted_x"] > cutoff)  | (median_df["Alexa488Adjusted_y"] > cutoff)]
        below_cut = median_df[~median_df.isin(above_cut)].dropna()
        print("Number of distinct names after filtering %d" % above_cut["Name"].nunique())

        above_cut = df_comb.merge(above_cut[["Name"]], on=["Name"])
        lbled = above_cut.groupby("Name").apply(lambda x :
            pd.Series({
                "p_coop": st.wilcox(
                        x["Alexa488Adjusted_y"].tolist(),
                        x["Alexa488Adjusted_x"].tolist(),
                        "greater"
                    ),
                "p_anti": st.wilcox(
                        x["Alexa488Adjusted_y"].tolist(),
                        x["Alexa488Adjusted_x"].tolist(),
                        "less"
                    ),
            })).reset_index()
        # too low statistical power, can't do fdr correction
        # lbled["p_coop"] = sm.fdrcorrection(lbled["p_coop"])[1]
        # lbled["p_anti"] = sm.fdrcorrection(lbled["p_anti"])[1]
        #print(lbled.groupby("p_coop").count())
        lbled["label"] = lbled.apply(
            lambda x:
                "cooperative" if x["p_coop"] < 0.061 else
                "anticoop" if x["p_anti"] < 0.061 else
                "additive",
                axis=1
            )
        df_lbled = median_df.merge(lbled, on="Name")
        # below_cut["label"] = "below_cutoff"
        df_lbled = df_lbled.append(below_cut, ignore_index=True)
        olist.append(df_lbled)
        arr.plot_classified_labels(df_lbled, path="%s_normalized.png"%ori, title="Cooperative plot, orientation %s" % ori, xlab=ch_x, ylab=ch_y)
        # write the labeled probes to csv
        df_wori = df_genomics[0][df_genomics[0]["ori"] == ori][["Name", "Sequence", "ets_pos", "runx_pos", "ori"]].drop_duplicates()
        #signif_label = df_lbled.loc[df_lbled["label"] != "below_cutoff",["Name","label","p_coop"]].drop_duplicates()
        lbled = df_lbled[["Name","label","p_coop"]].drop_duplicates()
        df_comb \
            .rename(columns={"Alexa488Adjusted_x":"ch1", "Alexa488Adjusted_y": "ch2"}) \
            .merge(lbled, on=["Name"]) \
            .merge(df_wori, on=["Name"]) \
            .to_csv("probes_labeled_%s.csv" % ori, index=False, float_format='%.4f')
    both_ori = olist[0].merge(olist[1], on=["Name"])

    print("Saving both orientation label")
    # this doesn't make sense, made just for the sake of having a combined plot
    both_ori["Alexa488Adjusted_x"] = (both_ori["Alexa488Adjusted_x_x"] + both_ori["Alexa488Adjusted_x_y"]) / 2
    both_ori["Alexa488Adjusted_y"] = (both_ori["Alexa488Adjusted_y_x"] + both_ori["Alexa488Adjusted_y_y"]) / 2
    both_ori["label"] = both_ori.apply(lambda x:
        "below_cutoff" if x["label_x"] == "below_cutoff" or x["label_y"] == "below_cutoff" else
        #"cooperative" if x["label_x"] == "cooperative" and x["label_y"] == "cooperative" else
        "anticoop" if x["p_coop_x"] == 1 and x["p_coop_y"] == 1 else
        "cooperative" if (x["p_coop_x"] < 0.061 and x["p_coop_y"] < 0.11) or (x["p_coop_y"] < 0.11 and x["p_coop_x"] < 0.61) else
        "additive",
        axis=1
    )
    both_ori[["Name","Alexa488Adjusted_x","Alexa488Adjusted_y","label"]].rename(columns={"Alexa488Adjusted_x":"ch1", "Alexa488Adjusted_y": "ch2"}).to_csv("labeled_w_medaff.csv", index=False, float_format='%.3f')
    both_ori_plt = both_ori[["Name","Alexa488Adjusted_x","Alexa488Adjusted_y","label"]]
    arr.plot_classified_labels(both_ori_plt, path="both_normalized.png", title="Cooperative plot, both orientations", xlab=ch_x, ylab=ch_y)

    print("Count per label:")
    print("er",olist[0][["label","Name"]].groupby("label").count())
    print("re",olist[1][["label","Name"]].groupby("label").count())
    print("both",both_ori_plt[["label","Name"]].groupby("label").count())

    name_info = df_genomics[0][["Name","Sequence","ets_pos","runx_pos","ori"]].drop_duplicates()
    all_labeled = name_info.merge(both_ori_plt[["Name","label"]], on="Name")
    all_labeled = all_labeled[all_labeled["label"] != "below_cutoff"]
    all_labeled["distance"] = abs(all_labeled["ets_pos"] - all_labeled["runx_pos"])
    all_labeled.to_csv("sequence_labeled_normalized.tsv", sep="\t", index=False, columns=["Name","Sequence","ets_pos","runx_pos","distance","ori","label"])

    print("Making boxplot of inconsistent label")
    # Make incosistency boxplot
    er, re = pd.read_csv("probes_labeled_er.csv"), pd.read_csv("probes_labeled_re.csv")
    indiv_df = pd.concat([er[["Name", "ori", "ch1"]], re[["Name", "ori", "ch1"]]]).rename(columns={'ch1':'affinity'})
    indiv_df.to_csv("indiv_df.tsv", sep="\t", index=False)
    two_df = pd.concat([er[["Name", "ori", "ch2"]], re[["Name", "ori", "ch2"]]]).rename(columns={'ch2':'affinity'})
    two_df.to_csv("two_df.tsv", sep="\t", index=False)
    lbldf = pd.merge(er[["Name","label"]].drop_duplicates(), re[["Name","label"]].drop_duplicates(), on="Name") \
            .rename(columns={'label_x':'er', 'label_y':'re'})
    lbldf.to_csv("lbl_df.tsv", sep="\t", index=False)
    arr.plot_ori_inconsistency(indiv_df, two_df, lbldf, log=True)

    print("Plotting additive p-value distribution")
    df_er = pd.read_csv("probes_labeled_er.csv")[["Name","p_coop","label"]] \
        .drop_duplicates(ignore_index=True)
    df_re = pd.read_csv("probes_labeled_re.csv")[["Name","p_coop","label"]] \
        .drop_duplicates(ignore_index=True)
    df = df_er.merge(df_re, on=["Name"], suffixes=("_er", "_re"))
    df_add_er = df[(df["label_er"] == "additive") & (df["label_re"] == "cooperative")][["p_coop_er"]].rename(columns={"p_coop_er":"p_coop"})
    df_add_er["type"] = "additive_er"
    df_add_re = df[(df["label_er"] == "cooperative") & (df["label_re"] == "additive")][["p_coop_re"]].rename(columns={"p_coop_re":"p_coop"})
    df_add_re["type"] = "additive_re"
    curdf = pd.concat([df_add_er, df_add_re])
    curdf["p_coop"] = curdf["p_coop"].astype(str)
    curdf.groupby('type')["p_coop"].value_counts().unstack(0).plot.barh()
    plt.title("p-val for additive in 1 orientation only")
    plt.savefig("additive_pval.png")
