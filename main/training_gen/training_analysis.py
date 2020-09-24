import os
os.chdir("../..")

import pandas as pd
import chip2probe.modeler.plotlib as plot
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

imads_ets_paths = ["input/site_models/imads_models/original/Ets1_20bp_GGAA.model", "input/site_models/imads_models/original/Ets1_20bp_GGAT.model"]
imads_ets_cores = ["GGAA", "GGAT"]
imads_ets_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_ets_paths, imads_ets_cores)]
imads_ets = iMADS(imads_ets_models, 0.19) # 0.2128 0.3061
#pickle.dump(imads_ets, open( "imads_ets_12.p", "wb" ) )

imads_runx_paths = ["input/site_models/imads_models/original/Runx1_20bp_GAGGT.model",
                "input/site_models/imads_models/original/Runx1_20bp_GCGGC.model",
                "input/site_models/imads_models/original/Runx1_20bp_GCGGG.model",
                "input/site_models/imads_models/original/Runx1_20bp_GTGGC.model",
                "input/site_models/imads_models/original/Runx1_20bp_GTGGG.model",
                "input/site_models/imads_models/original/Runx1_20bp_GTGGT.model"
                ]
imads_runx_cores = ["GAGGT", "GCGGC", "GCGGG", "GTGGC", "GTGGG", "GTGGT"]
imads_runx_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_runx_paths, imads_runx_cores)]
imads_runx = iMADS(imads_runx_models, 0.3) # 0.2128 0.3061
#pickle.dump(imads_runx, open( "imads_ets_20.p", "wb" ) )

def predict_strength(df, pred, tfname, flanklen=0):
    afflist = []
    for idx, row in df.iterrows():
        pos = row["%s_pos"%tfname]
        sites = pred.predict_sequence(row["fullseq"])
        aff = None
        for predsite in sites:
            predcpos = predsite["core_start"] - flanklen
            if pos > predcpos and pos < predcpos + predsite["core_width"]:
                aff = predsite['score']
                break
        afflist.append(aff)
    return afflist

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    train1_path = "%s/ch3_ch4/coop_ch3_vs_ch4/sequence_labeled.tsv" % basepath
    dftrain = pd.read_csv(train1_path, sep="\t").drop_duplicates(subset=["Name"], keep="first")
    fullseqs = pd.read_csv("%s/seqs_w_flanks.csv"% basepath)[["sequence","flank_left","flank_right"]].drop_duplicates(subset=["sequence"], keep="first")
    flanklen = len(fullseqs['flank_left'].iloc[0])

    dft = dftrain.merge(fullseqs, on="sequence")
    dft["fullseq"] = dft["flank_left"] + dft["sequence"] + dft["flank_right"]
    dft = dft.drop(["flank_left", "flank_right"], axis=1)

    dft["runx_score"] = predict_strength(dft, imads_runx, "runx", flanklen)
    dft["ets_score"] = predict_strength(dft, imads_ets, "ets", flanklen)

    dft = dft.drop(["fullseq"], axis=1)
    dft = dft[dft['ets_score'].notnull() & dft['runx_score'].notnull()]
    dft.to_csv("training.tsv",sep="\t",index=False)

    plot.plot_stacked_categories(dft, "distance")
    plot.plot_box_categories(dft, incols=["distance","runx_score", "ets_score"], alternative="smaller")
