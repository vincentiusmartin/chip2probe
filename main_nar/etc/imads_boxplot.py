import pandas as pd

import chip2probe.modeler.plotlib as plot

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

def make_imads_pred(df, imads):
    seqpred = []
    for idx, row in df.iterrows():
        rowdict = {"Name":row["Name"]}
        for tf in imads:
            preds = imads[tf].predict_sequence(row["Sequence"])
            for p in preds:
                if p["core_start"] == row["%s_start" % tf]:
                    rowdict["%s_score" % tf] = p["score"]
                    break
        if len(rowdict) == 3:
            rowdict["label"] = row["label"]
            seqpred.append(rowdict)
    return pd.DataFrame(seqpred)

pd.set_option("display.max_columns",None)
if __name__ == "__main__":
    # ========= iMADS =========
    basesp = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"

    imads_ets_paths = ["%s/input/sitemodels/imads_models/Ets1_w12_GGAA.model" % basesp, "%s/input/sitemodels/imads_models/Ets1_w12_GGAT.model" % basesp]
    imads_ets_cores = ["GGAA", "GGAT"]
    imads_ets_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_ets_paths, imads_ets_cores)]
    imads_ets = iMADS(imads_ets_models, 0.19) # 0.2128

    imads_runx_paths = ["%s/input/sitemodels/imads_models/Runx1_w20_GAGGT.model" % basesp,
                        "%s/input/sitemodels/imads_models/Runx1_w20_GCGGC.model" % basesp,
                        "%s/input/sitemodels/imads_models/Runx1_w20_GCGGG.model" % basesp,
                        "%s/input/sitemodels/imads_models/Runx1_w20_GCGGT.model" % basesp,
                        "%s/input/sitemodels/imads_models/Runx1_w20_GTGGC.model" % basesp,
                        "%s/input/sitemodels/imads_models/Runx1_w20_GTGGG.model" % basesp,
                        "%s/input/sitemodels/imads_models/Runx1_w20_GTGGT.model" % basesp]
    imads_runx_cores = ["GAGGT", "GCGGC", "GCGGG", "GCGGT", "GTGGC", "GTGGG", "GTGGT"]
    imads_runx_models = [iMADSModel(path, core, 20, [1, 2, 3]) for path, core in zip(imads_runx_paths, imads_runx_cores)]
    imads_runx = iMADS(imads_runx_models, 0.23) # 0.2128

    imads_dict = {"ets1":imads_ets, "runx1":imads_runx}

    # ==================

    # df = pd.read_csv("../output/Ets1Runx1/training/train_ets1_runx1.tsv", sep="\t")
    df = pd.read_csv("../output/Runx1Ets1/training/train_runx1_ets1.tsv", sep="\t")
    maintf, cooptf = "runx1", "ets1"
    dfpred = make_imads_pred(df, imads_dict)
    dfpred.to_csv("withimads.csv",index=False)


    dfpred.rename(columns={'%s_score' % maintf: '%s strength\n(main TF)' % maintf.capitalize(), '%s_score' % cooptf: '%s strength\n(cooperator TF)' % cooptf.capitalize()}, inplace=True)
    plot.plot_box_categories(dfpred, path="boxplot.pdf", incols=["%s strength\n(main TF)" % maintf.capitalize(), "%s strength\n(cooperator TF)" % cooptf.capitalize()], alternative="smaller", color=["#0343df","#75bbfd"])
