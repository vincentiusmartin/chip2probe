import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    ch12 = pd.read_csv("%s/ch1_ch2/coop_ch1_vs_ch2/labeled_w_medaff.csv" % basepath)
    ch34 = pd.read_csv("%s/ch3_ch4/coop_ch3_vs_ch4/labeled_w_medaff.csv" % basepath)

    ch12_in_ch34 = ch34[ch34["label"] == "cooperative"][["Name"]] \
        .merge(ch12, on="Name")

    ax = plt.axes()
    arr.plot_classified_labels(ch12, path="both.png",
                title="Cooperative plot, both orientations", col1="ch1", col2="ch2", axes=ax)

    x, y = np.log(ch12_in_ch34["ch1"].values), np.log(ch12_in_ch34["ch2"].values)
    ax.scatter(x, y, color="orange", s=3, label="coop_in_ch3_4")
    plt.savefig("bla.png")
