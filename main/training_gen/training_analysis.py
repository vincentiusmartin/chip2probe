import pandas as pd
import chip2probe.modeler.plotlib as plot

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    train1_path = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/labeled/ch1_vs_ch2/both_ori_labeled.tsv"
    df = pd.read_csv(train1_path, sep="\t").drop_duplicates(subset=["Name"], keep="first")
    plot.plot_stacked_categories(df, "distance")
    plot.plot_box_categories(df, incols=["distance"])
