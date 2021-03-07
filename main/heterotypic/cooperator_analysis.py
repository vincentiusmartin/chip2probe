import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import matplotlib.pyplot as plt
import numpy as np
import chip2probe.training_gen.arranalysis as arr

def analysis_overlap(df1, df2):
    print("Overlap in wtmt and ch1ch2")
    probe_overlap = df1.merge(df2, on=["Name","label"])
    print(probe_overlap["label"].value_counts())

    print("Only in ch1ch2")
    df1_only_names = set(df1["Name"]).difference(df2["Name"])
    df1_only_df = df1[df1["Name"].isin(df1_only_names)]
    print(df1_only_df["label"].value_counts())

    print("Only in wtmt")
    df2_only_names = set(df2["Name"]).difference(df1["Name"])
    df2_only_df = df2[df2["Name"].isin(df2_only_names)]
    print(df2_only_df["label"].value_counts())

if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    df1 = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/tables/labeled_w_medaff.csv" % basepath)
    df2 = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/tables/labeled_w_medaff.csv" % basepath)[["Name","label"]]
    df1 = df1[df1["label"] != "below_cutoff"]
    df2 = df2[df2["label"] != "below_cutoff"]
    analysis_overlap(df1, df2)

    df1_in_df2 = df2[df2["label"] == "cooperative"] \
        .merge(df1, on=["Name","label"])

    ax = plt.axes()
    arr.plot_classified_labels(df1, path="both.png",
                title="Cooperative plot, both orientations", col1="ch1", col2="ch2", axes=ax)
    x, y = np.log(df1_in_df2["ch1"].values), np.log(df1_in_df2["ch2"].values)
    ax.scatter(x, y, color="cyan", s=3, label="overlap coop ch1ch2 ch3ch4")
    ax.legend(loc="lower right")
    plt.savefig("bla.png")
