import pandas as pd
import pickle
import numpy as np

import chip2probe.util.stats_r as st
import statsmodels.stats.multitest as sm
from chip2probe.sitespredict.kompas import Kompas
import chip2probe.training_gen.arranalysis as arr

def get_single_site(seq, predictor):
    s = predictor.predict_sequence(seq)
    if len(s) == 1:
        return s[0]['core_mid']
    else:
        return None

def get_sites_pos(df, pred1, pred2, tf1, tf2, joincols=[], seqcol="Sequence"):
    """
    Get site position for each sequence

    Args:
        df: input data frame

    """
    seqdf = df[joincols+[seqcol]].drop_duplicates()
    seqdf["%s_pos"%tf1] = seqdf.apply(lambda x: get_single_site(x["Sequence"], pred1),axis=1)
    seqdf["%s_pos"%tf2] = seqdf.apply(lambda x: get_single_site(x["Sequence"], pred2),axis=1)
    seqdf = seqdf[~seqdf['%s_pos'%tf1].isna() & ~seqdf['%s_pos'%tf2].isna()]
    return seqdf

def read_probe_data(df, keyword, kompas_ets, kompas_runx):
    df_gen = df[df["Name"].str.contains(keyword,na=False)].sort_values(by=["Name"])
    df_gen[["Name","type","rep","ori"]] = df_gen["Name"].str.rsplit("_", n = 3, expand = True)
    df_gen["Sequence"] = df_gen["Sequence"].str[0:36]
    df_pos = get_sites_pos(df_gen[df_gen["type"] == "wt"], kompas_ets, kompas_runx, "ets", "runx", joincols=["Name","ori"])
    del df_pos["Sequence"]
    df_probe = df_pos.merge(df_gen, on=["Name","ori"])[["Name","Sequence","Alexa488Adjusted","type","ets_pos","runx_pos","ori","rep"]] \
                     .rename(columns={"Alexa488Adjusted":"affinity"})
    df_probe["ori"] = df_probe.apply(lambda x: "er" if x['ets_pos'] < x["runx_pos"] else "re", axis=1)
    return df_probe.sort_values(by=["Name","ori","rep"])

def assign_class(p, prevlbl):
    if prevlbl == "below_cutoff":
        return 'below_cutoff'
    elif prevlbl == 'anticoop' and p < 0.05:
        return 'anticoop'
    elif prevlbl == 'cooperative' and p < 0.05:
        return 'cooperative'
    else:
        return 'additive'

if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/coop_hetero-PBM_Ets_EtsRunx_v1"
    df = pd.read_csv("%s/Ets1_70.txt"%basepath,sep="\t") # Ets1_Runx1_70.txt
    cutoff = 651.4 # from plot_chamber_corr.py

    kompas_ets = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)
    kompas_runx = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Runx1_kmer_alignment.txt",
                    core_start = 12, core_end = 17, core_center = 14)

    pd.set_option("display.max_columns",None)

    # read the csv
    """
    keyword = "all_clean_seqs"
    print("Read probe data")
    df = read_probe_data(df, keyword, kompas_ets, kompas_runx)
    df.astype({'ets_pos': 'int32','runx_pos': 'int32'}).to_csv("df_wtmt_wpos.csv", index=False, float_format='%.3f')
    """
    # dflblpath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/df_labeled/wt_vs_mt/df_wtmt_wpos.csv"
    df = pd.read_csv("df_wtmt_wpos.csv")
    median_dict = df.groupby(["Name", "ori", "type"])["affinity"].median().to_dict()

    indivsum, twosites = arr.make_replicas_permutation(df)
    pickle.dump( indivsum, open( "indivsum.p", "wb" ) )
    pickle.dump( twosites, open( "twosites.p", "wb" ) )
    # indivsum = pickle.load( open( "indivsum.p", "rb" ) )
    # twosites = pickle.load( open( "twosites.p", "rb" ) )

    labeled_dict = {}
    for ori in ["er", "re"]:
        orilbls = []
        for k in indivsum[ori]:
            rowdict = {}
            if median_dict[(k,ori,'m3')] > cutoff or median_dict[(k,ori,'wt')] < cutoff:
                rowdict['label'], rowdict['p'] = "below_cutoff", 1
            else:
                rowdict['label'], rowdict['p'] =  arr.create_cooplbl(twosites[ori][k], indivsum[ori][k])
            rowdict['indiv_median'] = median_dict[(k,ori,'m1')] + median_dict[(k,ori,'m2')] - 2*median_dict[(k,ori,'m3')]
            rowdict['two_median'] = median_dict[(k,ori,'wt')] - median_dict[(k,ori,'m3')]
            rowdict['Name'] = k
            orilbls.append(rowdict)
        labeled_dict[ori] = pd.DataFrame(orilbls)
        print(labeled_dict[ori])
        labeled_dict[ori]['p'] = sm.fdrcorrection(labeled_dict[ori]['p'])[1]
        labeled_dict[ori]['label'] = labeled_dict[ori].apply(lambda row: assign_class(row['p'],row['label']),axis=1)
        print(ori, labeled_dict[ori]["label"].value_counts())
        arr.plot_classified_labels(labeled_dict[ori], col1="indiv_median", col2="two_median", log=True,
                        xlab="log(wt-m3)", ylab="log(m1-m3+m2-m3)", path="coop_log_%s.eps"%ori, title="Cooperative plot (in log), ori %s"%ori)
        arr.plot_classified_labels(labeled_dict[ori], col1="indiv_median", col2="two_median", log=False,
                        xlab="wt-m3", ylab="m1-m3+m2-m3", path="coop_%s.eps"%ori, title="Cooperative plot, ori %s"%ori)

    labeled_dict = pickle.load( open( "labeled.p", "rb" ) )
    df_er = labeled_dict['er'][labeled_dict['er']["label"] != "below_cutoff"]
    df_re = labeled_dict['re'][labeled_dict['re']["label"] != "below_cutoff"]
    labeled_joined = labeled_dict['er'].merge(labeled_dict['re'], on="Name", suffixes=('_er', '_re'))
    ori_match = labeled_joined[labeled_joined['label_er'] == labeled_joined['label_re']].rename(columns={'label_er':'label'})
    print("both", ori_match["label"].value_counts())
    ori_match['indiv_median'] = (ori_match['indiv_median_er'] + ori_match['indiv_median_re'])/2
    ori_match['two_median'] = (ori_match['two_median_er'] + ori_match['two_median_re'])/2
    arr.plot_classified_labels(ori_match, col1="indiv_median", col2="two_median", log=False,
                    xlab="wt-m3", ylab="m1-m3+m2-m3", path="coop_both.eps", title="Cooperative plot, both ori")
    arr.plot_classified_labels(ori_match, col1="indiv_median", col2="two_median", log=True,
                    xlab="wt-m3", ylab="m1-m3+m2-m3", path="coop_both_log.eps", title="Cooperative plot, both ori (log val)")

    #ori_match[["Name","label"]].to_csv("name_labeled.csv", index=False)
    labeled_per_ori = labeled_joined[["Name","label_er", "label_re"]]
    labeled_per_ori = labeled_per_ori[(labeled_per_ori["label_er"] != "below_cutoff") & (labeled_per_ori["label_re"] != "below_cutoff")] \
            .rename(columns={"label_er":"er","label_re":"re" })
    indivsum_ld = [{"Name":k, "ori": ori, "affinity":aff} for ori in indivsum for k in indivsum[ori] for aff in indivsum[ori][k]]
    indivsum_df = pd.DataFrame(indivsum_ld)
    twosites_ld = [{"Name":k, "ori": ori, "affinity":aff} for ori in twosites for k in twosites[ori] for aff in twosites[ori][k]]
    twosites_df = pd.DataFrame(twosites_ld)
    arr.plot_ori_incosistency2(indivsum_df, twosites_df, labeled_per_ori, log=True, fixed_ax=True)
