import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import numpy as np
import chip2probe.util.stats as st

if __name__=="__main__":
    pd.set_option('display.max_columns',None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/ch3_ch4/coop_ch3_vs_ch4/tables"
    seqmap = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/ch3_ch4/training_pwm.tsv",sep="\t").drop_duplicates()
    nm = "seq3143_all_clean_seqs_gap0_no20_0.4"
    #red/firebrick, mistyrose, blue, skyblue

    indiv = pd.read_csv("%s/indiv_df.tsv" % basepath, sep="\t")
    two = pd.read_csv("%s/two_df.tsv" % basepath, sep="\t")
    lbl = pd.read_csv("%s/lbl_df.tsv" % basepath, sep="\t")
    indiv["type"] = "one"
    two["type"] = "two"

    indiv = indiv.merge(lbl[["Name"]], on="Name")
    two = two.merge(lbl[["Name"]], on="Name")
    indiv["affinity"] = np.log(indiv["affinity"])
    two["affinity"] = np.log(two["affinity"])

    med_indiv = indiv.groupby("Name")[["affinity"]].median().reset_index()
    med_two = two.groupby("Name")[["affinity"]].median().reset_index()
    indiv_n= med_indiv[(med_indiv["affinity"] > 7.0) & (med_indiv["affinity"] < 8.5)]["Name"].unique().tolist()
    two_n = med_two[(med_two["affinity"] > 9.5)]["Name"].unique().tolist()
    print(set(indiv_n) & set(two_n))

    curdf = pd.concat([indiv, two])
    curdf = curdf[(curdf["Name"] == nm) & (curdf["ori"] == "er")]
    p_coop = st.wilcox(curdf[curdf["type"]=="two"]["affinity"].tolist(), curdf[curdf["type"]=="one"]["affinity"].tolist(), "greater")
    print("pc",p_coop)

    print("cd",curdf)
    print(seqmap[seqmap["Name"] == nm])

    arr.scatter_boxplot_col(curdf, "type", "affinity" , ["one","two"])
