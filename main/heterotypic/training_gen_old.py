import os
os.chdir("../..")

import pandas as pd
import chip2probe.modeler.plotlib as plot
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

import seaborn as sns
import matplotlib.pyplot as plt
from  matplotlib.ticker import FuncFormatter

from chip2probe.util import bio

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

pwm_ets = PWM("input/site_models/pwm/ets1.txt", log=True)
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
        if len(seq) < (corelen + 2*flanklen): # error
            res.append(-999)
            ori.append(-999)
        else:
            pred = pwm.predict_sequence(seq,zero_thres=False)[0]
            res.append(pred["score"])
            ori.append(pred["orientation"])
    return res, ori

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    train1_path = "%s/ch1_ch2/coop_ch1_vs_ch2/tables/sequence_labeled_normalized.tsv" % basepath
    dftrain = pd.read_csv(train1_path, sep="\t").drop_duplicates(subset=["Name"], keep="first")
    dftrain = dftrain[(dftrain["label"] == "cooperative") | (dftrain["label"] == "independent")]
    fullseqs = pd.read_csv("%s/seqs_w_flanks.csv"% basepath)[["sequence","flank_left","flank_right"]].drop_duplicates(subset=["sequence"], keep="first")
    flanklen = len(fullseqs['flank_left'].iloc[0])
    dft = dftrain.merge(fullseqs, on="sequence")
    dft["fullseq"] = dft["flank_left"] + dft["sequence"] + dft["flank_right"]
    dft = dft.drop(["flank_left", "flank_right"], axis=1)

    fullpos = pd.read_csv("output/heterotypic/EtsRunx_v1/seqpositions.csv")[["sequence", "ets_start", "runx_start"]]

    fullpos["ets_score"], fullpos["ets_ori"] = pwm_score(fullpos, pwm_ets, "ets_start", 4, 3)
    fullpos["runx_score"], fullpos["runx_ori"] = pwm_score(fullpos, pwm_runx, "runx_start", 5, 2)
    df = pd.read_csv("%s/df_wtmt_wpos.csv" % basepath)[["Name","Sequence"]].rename(columns={"Sequence":"sequence"}).drop_duplicates()
    df.merge(fullpos, on="sequence")[["Name","ets_score","runx_score"]].drop_duplicates().sort_values(by=["Name"]).to_csv("pwm_allseq.csv",index=False)

    dft = dft.merge(fullpos, on="sequence")
    dft = dft.drop(columns=["fullseq","ori"]).drop_duplicates()
    dft = dft[dft["distance"] != 4]
    # flip ets runx position
    for index, row in dft.iterrows():
        if row['ets_pos'] > row['runx_pos']:
            dft.at[index,'sequence'] = bio.revcompstr(row["sequence"])
            dft.at[index,'ets_start'] = len(row["sequence"]) - row["ets_pos"] - 4
            dft.at[index,'ets_pos'] = dft.at[index,'ets_start'] + 1
            dft.at[index,'ets_ori'] = 1 if row['ets_ori'] == -1 else 0
            dft.at[index,'runx_start'] = len(row["sequence"]) - row["runx_pos"] - 5
            dft.at[index,'runx_pos'] = dft.at[index,'runx_start'] + 2
            dft.at[index,'runx_ori'] = 1 if row['runx_ori'] == -1 else 0
        else:
            dft.at[index,'ets_ori'] = 0 if row['ets_ori'] == -1 else 1
            dft.at[index,'runx_ori'] = 0 if row['ets_ori'] == -1 else 1

    orimap = {0:"-",1:"+"}
    #dft["orientation"] = dft.apply(lambda x: "%s%s" % (str(x["ets_ori"]), str(x["runx_ori"])),axis=1)
    dft["orientation"] = dft.apply(lambda x: "%s/%s" % (orimap[int(x["ets_ori"])], orimap[int(x["runx_ori"])]),axis=1)
    #print(dft)
    dft.to_csv("training_pwm.tsv",sep="\t", index=False, float_format='%.3f')
    dft.rename(columns={'runx_score': 'Runx1 strength\n(cooperator TF)', 'ets_score': 'Ets1 strength\n(main TF)'}, inplace=True)
    #plot.plot_stacked_categories(dft, "distance", path="distance_bar.png", title="Distance distribution", ratio=True)
    plot.plot_stacked_categories(dft, "orientation", path="ori_bar.png", title="Relative sites orientation distribution", ratio=True)
    plot.plot_box_categories(dft, incols=["Runx1 strength\n(cooperator TF)", "Ets1 strength\n(main TF)"], alternative="smaller")

    # dft["runx_score"] = predict_strength(dft, pwm_runx, "runx", flanklen)
    # dft["ets_score"] = predict_strength(dft, pwm_ets, "ets", flanklen)
    # dft = dft.drop(["fullseq"], axis=1)
    # dft = dft[dft['ets_score'].notnull() & dft['runx_score'].notnull()]
    # dft.to_csv("training.tsv",sep="\t",index=False)
    # plot.plot_stacked_categories(dft, "distance")
    # plot.plot_box_categories(dft, incols=["distance","runx_score", "ets_score"], alternative="smaller")
