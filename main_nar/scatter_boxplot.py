import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import numpy as np
import chip2probe.util.stats as st


if __name__=="__main__":
    pd.set_option('display.max_columns',None)
    basepath = "output/Runx1Ets1/label_pr"
    # basepath = "output/Ets1Ets1/label_pr/"
    nm = "seq1790_all_clean_seqs_gap0_no20_0.4"
    ori = "er"

    indiv = pd.read_csv("%s/runx1_ets1_main.csv" % basepath)
    two = pd.read_csv("%s/runx1_ets1_main_cooperator.csv" % basepath)
    indiv["type"] = "one"
    two["type"] = "two"

    indiv["intensity"] = np.log(indiv["intensity"])
    two["intensity"] = np.log(two["intensity"])

    curdf = pd.concat([indiv, two])
    curdf = curdf[(curdf["Name"] == nm) & (curdf["ori"] == ori)]
    p_coop = st.wilcox(curdf[curdf["type"]=="two"]["intensity"].tolist(), curdf[curdf["type"]=="one"]["intensity"].tolist(), "greater")
    print("pc",p_coop)

    print("cd",curdf)

    arr.scatter_boxplot_col(curdf, "type", "intensity" , ["one","two"])
