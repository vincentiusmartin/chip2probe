import pandas as pd
from chip2probe.sitespredict.kompas import Kompas
import chip2probe.modeler.plotlib as plot
import numpy as np

def get_one_site_df(df, predictor, tfname, intcol="intensity"):
    seqdf = df[["Sequence"]].drop_duplicates()
    seqdf = seqdf[seqdf.apply(lambda x: len(predictor.predict_sequence(x["Sequence"])) == 1, axis=1)]
    outdf = df.merge(seqdf,on="Sequence").groupby(["Name"])["intensity"].median().reset_index().rename(columns={"intensity": "%s_intensity"%tfname})
    outdf["%s_intensity"%tfname] = np.log(outdf["%s_intensity"%tfname])
    return outdf

pd.set_option("display.max_columns",None)
kompas_ets = Kompas("../input/sitemodels/Ets1_kmer_alignment.txt", core_start = 11, core_end = 15, core_center = 12)
kompas_runx = Kompas("../input/sitemodels/Runx1_kmer_alignment.txt", core_start = 12, core_end = 17, core_center = 14)
if __name__ == "__main__":
    ets_raw = pd.read_csv("../input/probefiles/Ets1_only_pr_clean.csv")
    ets_mut = ets_raw[(ets_raw["type"] == "m1") | (ets_raw["type"] == "m2")]
    ets_df = get_one_site_df(ets_mut, kompas_ets, "ets")
    ets_df.to_csv("ets_one_site.csv",index=False)

    runx_raw = pd.read_csv("../input/probefiles/Runx1_only_pr_clean.csv")
    runx_mut = runx_raw[(runx_raw["type"] == "m1") | (runx_raw["type"] == "m2")]
    runx_df = get_one_site_df(runx_mut, kompas_runx, "runx")
    runx_df.to_csv("runx_one_site.csv",index=False)

    train_df = pd.read_csv("../output/Ets1Runx1/training/train_ets1_runx1.tsv", sep="\t")
    maintf, cooptf = "runx", "ets"
    color = ["#0343df","#75bbfd"] if maintf == "runx" else ["#b22222","#FFA07A"]

    intdf = ets_df.merge(runx_df, on="Name").merge(train_df[["Name","distance","orientation","label"]],on="Name")
    intdf.to_csv("intensity_ets_runx.csv",index=False)

    import sys
    sys.exit(0)

    intdf.rename(columns={'%s_intensity' % maintf: '%s intensity\n(main TF)' % maintf.capitalize(), '%s_intensity' % cooptf: '%s intensity\n(cooperator TF)' % cooptf.capitalize()}, inplace=True)
    plot.plot_box_categories(intdf, path="boxplot.pdf", incols=["%s intensity\n(main TF)" % maintf.capitalize(), "%s intensity\n(cooperator TF)" % cooptf.capitalize()], alternative="smaller", color=color)
