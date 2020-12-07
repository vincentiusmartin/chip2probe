import pandas as pd
import chip2probe.modeler.plotlib as plot
import matplotlib.pyplot as plt
import numpy as np

import chip2probe.training_gen.arranalysis as arr

def plot_sctr(df, col1, col2, path, title, xlab, ylab, coop_x, coop_y):
    plt.clf()
    ax = plt.axes()
    df['label_x'] = df['label_x'].replace({'cooperative': coop_x})
    coop_df = df[df["label_y"] == "cooperative"]
    arr.plot_classified_labels(df, col1=col1, col2=col2, labelcol="label_x",  title=title, xlab=xlab, ylab=ylab, labelnames=[coop_x,"independent","anticoop"],axes=ax)
    x, y = np.log(coop_df[col1].values), np.log(coop_df[col2].values)
    ax.scatter(x, y, color="blue", s=1, label=coop_y)
    ax.legend(loc="lower right")
    plt.savefig(path)

def combine_ch_intensity(df1, df2):
    one_df = df1[["Name","ch1","label"]].merge(df2[["Name","ch2","label"]], on="Name")
    two_df = df1[["Name","ch1","label"]].merge(df2[["Name","ch2","label"]], on="Name")
    plot_sctr(one_df, "ch1", "ch2", "corr_one.png", "Main TF only", "Ets1 only", "Runx1 only", "ER-coop", "RE-coop")
    plot_sctr(two_df, "ch1", "ch2", "corr_two.png", "Main TF + cooperator", "Ets1-Runx1", "Runx1-Ets1", "ER-coop", "RE-coop")


if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"

    df1_label = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/tables/both_ori_plt.csv" % basepath) \
            .rename(columns={"Alexa488Adjusted_x":"ch12_x","Alexa488Adjusted_y":"ch12_y"})[["Name","label"]]
    df2_label = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/tables/both_ori_plt.csv" % basepath) \
            .rename(columns={"Alexa488Adjusted_x":"ch34_x","Alexa488Adjusted_y":"ch34_y"})[["Name","label"]]

    df1_er = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/tables/probes_labeled_er.csv" % basepath)[["Name", "ch1", "ch2", "ori"]]
    df1_re = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/tables/probes_labeled_re.csv" % basepath)[["Name", "ch1", "ch2", "ori"]]
    df1 = pd.concat([df1_er, df1_re]).groupby(["Name","ori"]).median().reset_index().merge(df1_label, on="Name")

    df2_er = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/tables/probes_labeled_er.csv" % basepath)[["Name", "ch1", "ch2", "ori"]]
    df2_re = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/tables/probes_labeled_re.csv" % basepath)[["Name", "ch1", "ch2", "ori"]]
    df2 = pd.concat([df2_er, df2_re]).groupby(["Name","ori"]).median().reset_index().merge(df2_label, on="Name")

    arr.plot_classified_labels(df2, path="aaa.png", col1="ch1", col2="ch2", labelcol="label",  title="Ets_Runx1", xlab="Ch1", ylab="Ch2", labelnames=["cooperative","independent","anticoop"], )

    #combine_ch_intensity(df1, df2)


    """
    train_ch1 = pd.read_csv("%s/ch1_ch2/training_pwm.tsv" % basepath, sep="\t")
    train_ch2 = pd.read_csv("%s/ch3_ch4/training_pwm.tsv" % basepath, sep="\t")

    print("Chamber 1-2")
    print(df1.groupby("label")["ets_score"].count())

    print("Chamber 3-4")
    print(df2.groupby("label")["runx_score"].count())

    df2_additive_names = df2[df2["label"] == "independent"][["Name"]]
    df_both_add = df1[df1["label"] == "independent"].merge(df2_additive_names, on=["Name"])
    df2_coop = df2[df2["label"] == "cooperative"]
    df1_coop = df1[df1["label"] == "cooperative"]

    dft = pd.concat([df1_coop, df2_coop, df_both_add]).drop_duplicates()
    print(dft.groupby("label")["ets_score"].count())
    dft.to_csv("train_all.tsv",sep="\t",index=False)

    plot.plot_stacked_categories(dft, "distance", path="distance_bar.png", title="Distance distribution", ratio=True)
    plot.plot_stacked_categories(dft, "orientation", path="ori_bar.png", title="Relative sites orientation\ndistribution", ratio=True)
    plot.plot_box_categories(dft, incols=["runx_score","ets_score"], alternative="smaller")
    """
