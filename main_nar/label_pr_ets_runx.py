import chip2probe.training_gen.arranalysis as arr
import pandas as pd
import numpy as np
from chip2probe.sitespredict.kompas import Kompas

import chip2probe.util.stats_r as st

# =======

def get_one_site_seq(df, predictor):
    seqdf = df[["Sequence"]].drop_duplicates()
    seqdf["has_site"] = seqdf.apply(lambda x: len(predictor.predict_sequence(x["Sequence"])) == 1, axis=1)
    return seqdf

def get_med_m1m2(df, seqbound):
    d =  df[(df["type"] == "m1") | (df["type"] == "m2")]
    d = d.merge(seqbound, on="Sequence")
    d =  d[d["has_site"] == True]
    d = d.groupby(["Name","ori"])[["intensity"]].median().reset_index()
    return d

def get_normalization_coefs(df1, df2, predictor):
    """
    Using m1|m2 which are sequences with sites mutated
    """
    oneseqdf = get_one_site_seq(df1, predictor) # pd.read_csv("input/seq_w1_ets_site.csv") #
    oneseqdf["Sequence"] = oneseqdf["Sequence"].str[:36]
    med1, med2 = get_med_m1m2(df1, oneseqdf), get_med_m1m2(df2, oneseqdf)
    arr.plot_chamber_corr(med1, med2, median=False,
                           extrajoincols=["ori"], path="m1m2_log.png", log=False,
                           title="Correlation between chambers", xlab="Main TF only", ylab="Main TF + Cooperator TF", valcol="intensity")
    medcomb = med1.merge(med2, on=["Name","ori"])[["Name","ori", "intensity_x", "intensity_y"]]
    slope, intercept = np.polyfit(medcomb["intensity_x"], medcomb["intensity_y"], 1)

    return slope, intercept

def get_negcutoff(negdf, tf_str, slope, intercept, percentile):
    negdf = negdf[negdf["Name"].str.contains(tf_str)]
    negdf["intensity"]  = list(((negdf["intensity"] - intercept)/slope).values)
    return negdf[["intensity"]].quantile(percentile)

# =======

def assign_ori_pos(df, pred1, pred2, tf1, tf2, o1name, o2name, seqcol="Sequence"):
    """
    need to assign orientation based on which tf is the main TF
    We put the main TF near the free end
    """
    dfgen = df[df["type"] == "wt"]
    del dfgen["ori"]
    seqdf = dfgen[[seqcol]].drop_duplicates()

    tfpred = [[tf1, pred1], [tf2, pred2]]
    for tf, pred in tfpred:
        poslist = []
        startlist = []
        for index,row in seqdf.iterrows():
            s = pred.predict_sequence(row["Sequence"])
            if len(s) == 1:
                poslist.append(s[0]['core_mid'])
                startlist.append(s[0]['core_start'])
            else:
                poslist.append(-999)
                startlist.append(-999)
        seqdf["%s_pos"%tf] = poslist
        seqdf["%s_start"%tf] = startlist

    seqdf["ori"] = seqdf.apply(lambda x: o1name if x["%s_pos"%tf1] < x["%s_pos"%tf2] else o2name,axis=1)
    seqdf = seqdf[~seqdf['%s_pos'%tf1].isna() & ~seqdf['%s_pos'%tf2].isna()]
    probedf = seqdf.merge(dfgen, on="Sequence")[["Name","Sequence","intensity","%s_pos" % tf1, "%s_start" % tf1, "%s_pos" % tf2, "%s_start" % tf2, "ori","rep"]]
    return probedf


# label the probes for ets1 ets1
kompas_ets = Kompas("input/sitemodels/Ets1_kmer_alignment.txt", core_start = 11, core_end = 15, core_center = 12)
kompas_runx = Kompas("input/sitemodels/Runx1_kmer_alignment.txt", core_start = 12, core_end = 17, core_center = 14)
if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    pd.set_option("mode.chained_assignment",None)
    neg_percentile = 0.75
    p_default = 0.015
    p_ambiguous = 0.06

    # For plotting
    # m for main TF, mc for main and cooperator TF
    ch_x =  "Ets1 only"
    ch_y = "Ets1 + Runx1"
    both_title = "Cooperative vs independent binding of Ets1-Runx1"
    maintf = "ets1"
    cooptf = "runx1"
    df_m, neg_m = pd.read_csv("input/probefiles/Ets1_only_pr_clean.csv"), pd.read_csv("input/probefiles/Ets1_only_neg_clean.csv")
    df_mc, neg_mc = pd.read_csv("input/probefiles/Ets1_Runx1_pr_clean.csv"), pd.read_csv("input/probefiles/Ets1_Runx1_neg_clean.csv")
    baseoutpath = "output/Ets1Runx1/label_pr/"

    # ch_x =  "Runx1 only"
    # ch_y = "Runx1 + Ets1"
    # both_title = "Cooperative vs independent binding of Runx1-Ets1"
    # maintf = "runx1"
    # cooptf = "ets1"
    # df_m, neg_m = pd.read_csv("input/probefiles/Runx1_only_pr_clean.csv"), pd.read_csv("input/probefiles/Runx1_only_neg_clean.csv")
    # df_mc, neg_mc = pd.read_csv("input/probefiles/Runx1_Ets1_pr_clean.csv"), pd.read_csv("input/probefiles/Runx1_Ets1_neg_clean.csv")
    # baseoutpath = "output/Runx1Ets1/label_pr/"

    print("Number of raw probes distinct names %d" % df_m["Name"].nunique())


    # if main tf ets1, we say cooperative everytime it's cooperative on re orientation
    label_er_re = ["independent", "cooperative"] if maintf == "ets1" else ["cooperative", "independent"]
    #ori_er_re = "re" if maintf == "ets1" else "er"
    main_predictor = kompas_ets if maintf == "ets1" else kompas_runx

    slope,intercept = get_normalization_coefs(df_m, df_mc, main_predictor)
    print("Normalization using y = %.2fx + %.2f" % (slope,intercept))

    cutoff = float(get_negcutoff(neg_m, cooptf, slope, intercept, neg_percentile))
    print("Negative control cutoff %.3f" % cutoff)

    df_m = assign_ori_pos(df_m, kompas_ets, kompas_runx, "ets1", "runx1" , "er", "re") # pd.read_csv("%s_%s_main.csv" % (maintf, cooptf)) #a
    df_mc = assign_ori_pos(df_mc, kompas_ets, kompas_runx, "ets1", "runx1", "er", "re") # pd.read_csv("%s_%s_main_cooperator.csv" % (maintf, cooptf))
    df_mc["intensity"] = (df_mc["intensity"] - intercept)/slope
    df_m.to_csv("%s/%s_%s_main.csv" % (baseoutpath, maintf, cooptf),index=False)
    df_mc.to_csv("%s/%s_%s_main_cooperator.csv" % (baseoutpath, maintf, cooptf),index=False)

    olist = []
    for ori in ["er","re"]:
        print("Orientation %s" % ori)
        cur_m = df_m[df_m["ori"] == ori][["Name","rep","intensity"]]
        cur_mc = df_mc[df_mc["ori"] == ori][["Name","rep","intensity"]]
        comb = cur_m.merge(cur_mc, on=["Name","rep"])

        print("Number of distinct names before filtering %d" % comb["Name"].nunique())
        median_df = comb.groupby(["Name"]).median().reset_index()
        above_cut = median_df[(median_df["intensity_x"] > cutoff)  | (median_df["intensity_y"] > cutoff)]
        below_cut = median_df[~median_df.isin(above_cut)].dropna()
        print("Number of distinct names after filtering %d" % above_cut["Name"].nunique())

        above_cut = comb.merge(above_cut[["Name"]], on=["Name"])
        lbled = above_cut.groupby("Name").apply(lambda x :
            pd.Series({
                "p_coop": st.wilcox(x["intensity_y"].tolist(), x["intensity_x"].tolist(), "greater"),
                "p_anti": st.wilcox(x["intensity_y"].tolist(), x["intensity_x"].tolist(), "less"),
            })
        ).reset_index()

        lbled["label"] = lbled.apply(
            lambda x:
                "cooperative" if x["p_coop"] < p_default else
                "ambiguous" if x["p_coop"] < p_ambiguous else
                "anticooperative" if x["p_anti"] < p_default else
                "independent",
                axis=1
        )

        df_lbled = median_df.merge(lbled, on="Name")
        below_cut["label"] = "below_cutoff"
        df_lbled = df_lbled.append(below_cut, ignore_index=True)
        olist.append(df_lbled)

    both_ori = olist[0].merge(olist[1], on=["Name"], suffixes=("_er", "_re"))
    print("Number of distinct names, above cutoff, after orientation joining" % both_ori[both_ori["label"] != "below_cutoff"]["Name"].nunique())

    # labeling the probe using both orientations
    oricoop = "er" if cooptf == "ets1" else "re"
    both_ori["intensity_x"] = both_ori["intensity_x_%s" % oricoop] #(both_ori["intensity_x_er"] + both_ori["intensity_x_re"]) / 2
    both_ori["intensity_y"] = both_ori["intensity_y_%s" % oricoop]  #(both_ori["intensity_y_er"] + both_ori["intensity_y_re"]) / 2
    both_ori["label"] = both_ori.apply(lambda x:
        "below_cutoff" if x["label_er"] == "below_cutoff" or x["label_re"] == "below_cutoff" else
        "anticooperative" if x["p_coop_er"] == 1 and x["p_coop_re"] == 1 else
        "cooperative" if (x["label_er"] == "cooperative" and x["label_re"] == "cooperative") or
                         (x["label_er"] == "cooperative" and x["label_re"] == "ambiguous") or
                         (x["label_er"] == "ambiguous" and x["label_re"] == "cooperative") or
                         (x["label_er"] == label_er_re[0] and x["label_re"] == label_er_re[1]) else
        "independent" if (x["label_er"] == "independent" and x["label_re"] == "independent") or
                    (x["label_er"] == "independent" and x["label_re"] == "ambiguous") or
                    (x["label_er"] == "ambiguous" and x["label_re"] == "independent") else
        "ambiguous",
        axis=1
    )
    # filter out anticoop since we only have a few
    both_ori = both_ori[both_ori["label"] != "anticooperative"]

    both_ori_plt = both_ori[["Name","intensity_x","intensity_y","label"]]
    both_ori_plt.to_csv("%s/both_ori_plt_%s_%s.csv" % (baseoutpath,maintf,cooptf),index=False)
    arr.plot_classified_labels(both_ori_plt, path="%s/both_normalized_%s_%s.png" % (baseoutpath,maintf,cooptf), col1="intensity_x", col2="intensity_y",
            title=both_title, xlab=ch_x, ylab=ch_y, plotnonsignif=False, labelnames=["cooperative","independent","anticoop"])

    print("Count per label",both_ori_plt[["label","Name"]].groupby("label").count())

    # we use sequence where ets1 is on the left for simplicity
    name_info = df_m[df_m["ori"] == "er"][["Name","Sequence","%s_pos" % maintf,"%s_start" % maintf, "%s_pos" % cooptf, "%s_start" % cooptf, "ori"]].drop_duplicates()
    all_labeled = name_info.merge(both_ori_plt[["Name","label"]], on="Name")
    all_labeled = all_labeled[all_labeled["label"] != "below_cutoff"]
    all_labeled["distance"] = abs(all_labeled["%s_pos" % maintf] - all_labeled["%s_pos" % cooptf])
    all_labeled.to_csv("%s/seqbled_%s_%s.tsv" % (baseoutpath,maintf,cooptf), sep="\t", index=False, columns=["Name", "Sequence", "%s_pos" % maintf, "%s_start" % maintf, "%s_pos" % cooptf, "%s_start" % cooptf, "distance", "ori", "label"])
