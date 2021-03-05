import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import pickle
import matplotlib.pyplot as plt
import numpy as np

def process_arrfile(df):
    df[["Name","type"]] = df["Name"].str.rsplit('_',n=1,expand=True)
    df = df.melt(id_vars=["Name", "type", "Sequence"], var_name="ori", value_name="affinity")
    df = df[~df["ori"].str.contains("Median")]
    df[["ori","rep"]] = df["ori"].str.split("_",expand=True)
    return df

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)

    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/homotypic/training"

    selected = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/array_design_files/Coop2Ets_validation/custom_probes_selected.csv")
    wtseqs = selected[selected["comment"] == "wt"][["sequence","wtlabel"]].drop_duplicates().rename({"sequence":"Sequence"}, axis=1)
    print(wtseqs["wtlabel"].value_counts())
    wtcoop = wtseqs[wtseqs['wtlabel'] == 1]
    wtadd = wtseqs[wtseqs['wtlabel'] == 0]
    snm = pd.read_csv("%s/seqnamemap.csv" % basepath)
    wtnames_coop = wtcoop.merge(snm, on="Sequence")[["Name"]].drop_duplicates()
    wtnames_add = wtadd.merge(snm, on="Sequence")[["Name"]].drop_duplicates()
    print(wtcoop)

    # for ori in ['o1','o2']:
    #     lo = pd.read_csv("%s/lbled_%s.csv" % (basepath,ori))
    #     ax = plt.axes()
    #     arr.plot_classified_labels(lo, col1="indiv_median", col2="two_median", log=True,
    #                     xlab="log(m1-m3+m2-m3)", ylab="log(wt-m3)", path="labeled_log_%s.png" % ori,
    #                     title="Cooperative plot, %s" % ori, axes=ax)
    #     print("Count %s" % ori,lo["label"].value_counts())
    #     wtcoop_aff = lo.merge(wtnames_coop)
    #     wtadd_aff = lo.merge(wtnames_add)
    #     ax.scatter(np.log(wtcoop_aff["indiv_median"]), np.log(wtcoop_aff["two_median"]), color="cyan", s=3, label="coop_in_cust")
    #     ax.scatter(np.log(wtadd_aff["indiv_median"]), np.log(wtadd_aff["two_median"]), color="green", s=3, label="add_in_cust")
    #     ax.legend(loc="lower right")
    #     plt.savefig("cust_in_ori_%s.png" % ori)
    #     plt.clf()
