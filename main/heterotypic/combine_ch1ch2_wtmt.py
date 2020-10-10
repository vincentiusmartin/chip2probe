
import pandas as pd
import math
import chip2probe.training_gen.arranalysis as arr
import matplotlib.pyplot as plt
import numpy as np

def label_probe(p1, p2, pcut1, pcut2):
    if math.isnan(p1) or math.isnan(p2):
        return "below_cutoff"
    elif p1 == 1 and p2 == 1:
        return "anticoop"
    elif (p1 < pcut1 and p2 < pcut2) or (p1 < pcut2 and p2 < pcut1):
        return "cooperative"
    else:
        return "additive"

def get_med_df(ch_er, ch_re):
    ch_med = ch_er[["Name","ch1","ch2"]].groupby(["Name"]).median().reset_index().merge(
        ch_re[["Name","ch1","ch2"]].groupby(["Name"]).median().reset_index(), on=["Name"]
    )
    ch_med["ch1"] = (ch_med["ch1_x"] + ch_med["ch1_y"]) / 2
    ch_med["ch2"] = (ch_med["ch2_x"] + ch_med["ch2_y"]) / 2
    ch_med = ch_med[["Name", "ch1", "ch2"]]
    return ch_med


def read_tbls(df_o1, df_o2, pcut1, pcut2):
    ch_er = df_o1[["Name","p_coop","ori"]].drop_duplicates()
    ch_re = df_o2[["Name","p_coop","ori"]].drop_duplicates()
    ch_med = get_med_df(df_o1[["Name","ch1","ch2"]].drop_duplicates(), df_o2[["Name","ch1","ch2"]].drop_duplicates())
    ch = ch_er.merge(ch_re, on="Name", suffixes=("_er","_re"))
    ch["label"] = ch.apply(lambda x : label_probe(x['p_coop_er'], x["p_coop_re"], pcut1, pcut2), axis=1)
    ch["label_er"] = ch.apply(lambda x: label_one_ori(x["p_coop_er"], pcut1), axis=1)
    ch["label_re"] = ch.apply(lambda x: label_one_ori(x["p_coop_re"], pcut1), axis=1)
    return ch.merge(ch_med, on="Name")

def get_diff(df1, df2, label):
    df1_spec = df1[df1["label"] == label]
    df2_spec = df2[df2["label"] == label]
    df1_only_names = set(df1_spec["Name"]) - set(df2_spec["Name"])
    print(len(df1_only_names))
    df1_only_df = df1[df1["Name"].isin(df1_only_names)]
    return df1_only_df

def label_one_ori(p, pcut):
    if math.isnan(p):
        return "below_cutoff"
    elif p == 1:
        return "anticoop"
    elif p < pcut:
        return "cooperative"
    else:
        return "additive"

if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    pd.set_option("display.max_columns",None)
    seqpwm = pd.read_csv("%s/pwm_allseq.csv"%basepath)

    ch12_er = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/tables/probes_labeled_er.csv" % basepath, sep=",")
    ch12_re = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/tables/probes_labeled_re.csv" % basepath, sep=",")
    ch12 = read_tbls(ch12_er, ch12_re, 0.061, 0.11)
    ch12_indiv = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/tables/indiv_df.tsv" % basepath,sep="\t")
    ch12_two = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/tables/two_df.tsv" % basepath,sep="\t")

    ch34_er = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/tables/probes_labeled_er.csv" % basepath, sep=",")
    ch34_re = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/tables/probes_labeled_re.csv" % basepath, sep=",")
    ch34 = read_tbls(ch34_er, ch34_re, 0.015, 0.015)
    ch34_indiv = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/tables/indiv_df.tsv" % basepath,sep="\t")
    ch34_two = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/tables/two_df.tsv" % basepath,sep="\t")

    print(ch12)
    print("ch1ch2")
    print(ch12["label"].value_counts())

    print("ch3ch4")
    print(ch34["label"].value_counts())

    # plot
    ax = plt.axes()
    arr.plot_classified_labels(ch12,
                title="Cooperative plot, both orientations", col1="ch1", col2="ch2", axes=ax)
    coop34_in12 = ch34[ch34["label"] == "cooperative"][["Name"]] \
        .merge(ch12, on=["Name"])
    x, y = np.log(coop34_in12["ch1"].values), np.log(coop34_in12["ch2"].values)
    ax.scatter(x, y, color="cyan", s=3, label="overlap coop ch1ch2 ch3ch4")
    ax.legend(loc="lower right")
    plt.savefig("overlap.png")

    print("Overlap in ch1ch2 and ch3ch4")
    probe_overlap = ch12.merge(ch34, on="Name")
    match_lbl = probe_overlap[probe_overlap["label_x"] == probe_overlap["label_y"]]
    print(probe_overlap.shape[0])
    print(match_lbl["label_x"].value_counts())

    print("Only in ch1ch2")
    coop12 = set(get_diff(ch12, ch34, "cooperative")["Name"])
    ch34_coop12 = ch34[ch34["Name"].isin(coop12)]
    lbl34_12 = ch34_coop12[["Name","label_er", "label_re"]].rename(columns={'label_er':'er', 'label_re':'re'}).drop_duplicates()
    print(ch34_coop12.groupby(['label_er','label_re']).size().reset_index().rename(columns={0:'count'}))
    # arr.plot_ori_inconsistency(ch34_indiv[ch34_indiv["Name"].isin(coop12)],
    #                            ch34_two[ch34_two["Name"].isin(coop12)], lbl34_12, log=True)

    print("Only in ch3ch4")
    coop34 = set(get_diff(ch34, ch12, "cooperative")["Name"])
    ch12_coop34 = ch12[ch12["Name"].isin(coop34)]
    lbl12_34 = ch12_coop34[["Name","label_er", "label_re"]].rename(columns={'label_er':'er', 'label_re':'re'}).drop_duplicates()
    print(ch12_coop34.groupby(['label_er','label_re']).size().reset_index().rename(columns={0:'count'}))
    # arr.plot_ori_inconsistency(ch12_indiv[ch12_indiv["Name"].isin(coop34)],
    #                             ch12_two[ch12_two["Name"].isin(coop34)], lbl12_34, log=True)

    plt.clf()
    dfplot = ch34_coop12[(ch34_coop12["label_er"] == "cooperative") & (ch34_coop12["label_re"] == "below_cutoff")][["Name"]].merge(seqpwm, on="Name")
    dfplot.boxplot(column=["ets_score","runx_score"])
    plt.title("PWM score--cooperative_er below_cutoff_re")
    plt.savefig("box.png")
