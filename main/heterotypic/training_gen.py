import os
os.chdir("../..")

import pandas as pd
import chip2probe.modeler.plotlib as plot
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

from chip2probe.sitespredict.pwm import PWM

# imads_ets_paths = ["input/site_models/imads_models/original/Ets1_20bp_GGAA.model", "input/site_models/imads_models/original/Ets1_20bp_GGAT.model"]
# imads_ets_cores = ["GGAA", "GGAT"]
# imads_ets_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_ets_paths, imads_ets_cores)]
# imads_ets = iMADS(imads_ets_models, 0.19) # 0.2128 0.3061
#
# imads_runx_paths = ["input/site_models/imads_models/original/Runx1_20bp_GAGGT.model",
#                 "input/site_models/imads_models/original/Runx1_20bp_GCGGC.model",
#                 "input/site_models/imads_models/original/Runx1_20bp_GCGGG.model",
#                 "input/site_models/imads_models/original/Runx1_20bp_GTGGC.model",
#                 "input/site_models/imads_models/original/Runx1_20bp_GTGGG.model",
#                 "input/site_models/imads_models/original/Runx1_20bp_GTGGT.model"
#                 ]
# imads_runx_cores = ["GAGGT", "GCGGC", "GCGGG", "GTGGC", "GTGGG", "GTGGT"]
# imads_runx_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_runx_paths, imads_runx_cores)]
# imads_runx = iMADS(imads_runx_models, 0.3) # 0.2128 0.3061

pwm_ets = PWM("input/site_models/pwm/ets1.txt", 6, 16, log=True)
pwm_runx = PWM("input/site_models/pwm/runx1.txt", 8, 17, log=True)

def predict_strength(df, pred, tfname, flanklen=0):
    afflist = []
    for idx, row in df.iterrows():
        pos = row["%s_pos"%tfname]
        sites = pred.predict_sequence(row["fullseq"])
        aff = None
        for predsite in sites:
            predcpos = predsite["site_start"] - flanklen
            if pos > predcpos and pos < predcpos + predsite["site_width"]:
                aff = predsite['score']
                break
        afflist.append(aff)
    return afflist

def pwm_score(dfpos, pwm, startcol, corelen, flanklen):
    core_start = pwm.length // 2 - corelen // 2
    core_end = core_start + corelen
    res = []
    ori = []
    for idx, row in dfpos.iterrows():
        seq = row["sequence"][row[startcol]-flanklen:row[startcol]+corelen+flanklen]
        #print(row["sequence"], seq, startcol, row[startcol])
        pred = pwm.predict_sequence(seq,zero_thres=False)[0]
        res.append(pred["score"])
        ori.append(pred["orientation"])
    return res, ori

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    train1_path = "%s/ch1_ch2/coop_ch1_vs_ch2/sequence_labeled_p06_p1.tsv" % basepath
    dftrain = pd.read_csv(train1_path, sep="\t").drop_duplicates(subset=["Name"], keep="first")
    fullseqs = pd.read_csv("%s/seqs_w_flanks.csv"% basepath)[["sequence","flank_left","flank_right"]].drop_duplicates(subset=["sequence"], keep="first")
    flanklen = len(fullseqs['flank_left'].iloc[0])
    dft = dftrain.merge(fullseqs, on="sequence")
    dft["fullseq"] = dft["flank_left"] + dft["sequence"] + dft["flank_right"]
    dft = dft.drop(["flank_left", "flank_right"], axis=1)

    fullpos = pd.read_csv("output/heterotypic/EtsRunx_v1/seqpositions.csv")[["sequence", "ets_start", "runx_start"]] \
        .merge(dft[["sequence"]], on="sequence")

    fullpos["ets_score"], fullpos["ets_ori"] = pwm_score(fullpos, pwm_ets, "ets_start", 4, 3)
    fullpos["runx_score"], fullpos["runx_ori"] = pwm_score(fullpos, pwm_runx, "runx_start", 5, 2)

    dft = dft.merge(fullpos, on="sequence")
    dft = dft.drop(columns=["fullseq","ori"]).drop_duplicates()
    dft.to_csv("training_pwm.tsv",sep="\t", index=False, float_format='%.3f')
    plot.plot_stacked_categories(dft, "distance")
    plot.plot_box_categories(dft, incols=["distance","runx_score", "ets_score"], alternative="smaller")

    """
    dft["runx_score"] = predict_strength(dft, pwm_runx, "runx", flanklen)
    dft["ets_score"] = predict_strength(dft, pwm_ets, "ets", flanklen)
    dft = dft.drop(["fullseq"], axis=1)
    dft = dft[dft['ets_score'].notnull() & dft['runx_score'].notnull()]
    dft.to_csv("training.tsv",sep="\t",index=False)
    plot.plot_stacked_categories(dft, "distance")
    plot.plot_box_categories(dft, incols=["distance","runx_score", "ets_score"], alternative="smaller")
    """
