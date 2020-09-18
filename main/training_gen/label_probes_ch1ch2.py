import pandas as pd
import chip2probe.util.stats_r as st

from chip2probe.sitespredict.kompas import Kompas
import chip2probe.training_gen.arranalysis as arr


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
    dflist = [pd.read_csv("%s/Ets1_70.txt"%basepath,sep="\t"),
            pd.read_csv("%s/Ets1_Runx1_70.txt"%basepath,sep="\t")]
    cutoff = 651.4 # from plot_chamber_corr.py

    kompas_ets = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)
    kompas_runx = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Runx1_kmer_alignment.txt",
                    core_start = 12, core_end = 17, core_center = 14)

    # ------- labeling sequences, for now only take wt -------
    keyword = "all_clean_seqs"
    df_genomics = []

    for i in range(len(dflist)):
        df = dflist[i]
        df_gen = read_probe_data(df)
        if i == 1:
          df_gen["Alexa488Adjusted"] = (df_gen["Alexa488Adjusted"] - 107.74) / 0.82 # normalize
        df_gen.to_csv("df_wt_ch%d.csv" % (i+1), index=False)
        df_gen["Sequence"] = df_gen["Sequence"].str[0:36]
        df_genomics.append(df_gen)
    df_genomics = [pd.read_csv("df_wt_ch1.csv").convert_dtypes(), pd.read_csv("df_wt_ch2.csv").convert_dtypes()]

    # we filter each group if median intensity is lower than cutoff
    olist = []

    for ori in ["er","re"]:
        print("checking orientation %s" % ori)
        df1 = df_genomics[0][df_genomics[0]["ori"] == ori][["Name","rep","Alexa488Adjusted",]]
        df2 = df_genomics[1][df_genomics[1]["ori"] == ori][["Name","rep","Alexa488Adjusted"]]
        df_comb = df1.merge(df2, on=["Name","rep"])
        print("Number of distinct names before filtering %d" % df_comb["Name"].nunique())
        median_df = df_comb.groupby(["Name"]).median().reset_index()
        above_cut = median_df[(median_df["Alexa488Adjusted_x"] > cutoff) & (median_df["Alexa488Adjusted_y"] > cutoff)]
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
        lbled["label"] = lbled.apply(
            lambda x:
                "cooperative" if x["p_coop"] < 0.06 else
                "anticoop" if x["p_anti"] < 0.06 else
                "additive",
                axis=1
            )
        df_lbled = median_df.merge(lbled, on="Name")
        #df_lbled = df_lbled.drop(['p_coop', 'p_anti'], axis=1)
        below_cut["label"] = "below_cutoff"
        df_lbled = df_lbled.append(below_cut, ignore_index=True)
        olist.append(df_lbled)
        arr.plot_classified_labels(df_lbled, filepath="%s_normalized.png"%ori, title="Cooperative plot, orientation %s" % ori)
        # write the labeled probes to csv
        df_wori = df_genomics[0][df_genomics[0]["ori"] == ori][["Name", "Sequence", "ets_pos", "runx_pos", "ori"]].drop_duplicates()
        signif_label = df_lbled.loc[df_lbled["label"] != "below_cutoff",["Name","label"]].drop_duplicates()
        df_comb \
            .rename(columns={"Alexa488Adjusted_x":"ch1", "Alexa488Adjusted_y": "ch2"}) \
            .merge(signif_label, on=["Name"]) \
            .merge(df_wori, on=["Name"]) \
            .to_csv("probes_labeled_%s.csv" % ori, index=False, float_format='%.4f')
    both_ori = olist[0].merge(olist[1], on=["Name"])

    # this doesn't make sense, made just for the sake of having a combined plot
    both_ori["Alexa488Adjusted_x"] = (both_ori["Alexa488Adjusted_x_x"] + both_ori["Alexa488Adjusted_x_y"]) / 2
    both_ori["Alexa488Adjusted_y"] = (both_ori["Alexa488Adjusted_y_x"] + both_ori["Alexa488Adjusted_y_y"]) / 2
    both_ori["label"] = both_ori.apply(lambda x:
        "below_cutoff" if x["label_x"] == "below_cutoff" or x["label_y"] == "below_cutoff" else
        "cooperative" if x["label_x"] == "cooperative" and x["label_y"] == "cooperative" else
        "anticoop" if x["label_x"] == "anticoop" and x["label_y"] == "anticoop" else
        "additive",
        axis=1
    )
    both_ori = both_ori[["Name","Alexa488Adjusted_x","Alexa488Adjusted_y","label"]]
    arr.plot_classified_labels(both_ori, filepath="both_normalized.png", title="Cooperative plot, both orientations")

    print(olist[0][["label","Name"]].groupby("label").count())
    print(olist[1][["label","Name"]].groupby("label").count())
    print(both_ori[["label","Name"]].groupby("label").count())

    name_info = df_genomics[0][["Name","Sequence","ets_pos","runx_pos","ori"]].drop_duplicates()
    all_labeled = name_info.merge(both_ori[["Name","label"]], on="Name")
    all_labeled = all_labeled[all_labeled["label"] != "below_cutoff"]
    all_labeled["distance"] = abs(all_labeled["ets_pos"] - all_labeled["runx_pos"])
    all_labeled.to_csv("sequence_labeled_normalized.tsv", sep="\t", index=False, columns=["Name","Sequence","ets_pos","runx_pos","distance","ori","label"])
