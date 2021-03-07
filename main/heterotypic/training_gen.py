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

pwm_ets = PWM("input/site_models/pwm/ets1.txt", log=True, reverse=False)
pwm_runx = PWM("input/site_models/pwm/runx1.txt", 8, 17, log=True, reverse=True)

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
    res = []
    ori = []
    core = []
    for idx, row in dfpos.iterrows():
        seq = row["sequence"][row[startcol]-flanklen:row[startcol]+corelen+flanklen]
        core.append(row['sequence'][row[startcol]:row[startcol]+corelen])
        if len(seq) < (corelen + 2*flanklen): # error
            res.append(-999)
            ori.append(-999)
        else:
            pred = pwm.predict_sequence(seq,zero_thres=False)[0]
            res.append(pred["score"])
            cur_ori = 0 if pred["orientation"] == -1 else 1
            ori.append(cur_ori)
    return res, ori, core


if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    train1_path = "%s/ch1_ch2/coop_ch1_vs_ch2/tables/sequence_labeled_normalized.tsv" % basepath
    dftrain = pd.read_csv(train1_path, sep="\t").drop_duplicates(subset=["Name"], keep="first")
    dftrain = dftrain[(dftrain["label"] == "cooperative") | (dftrain["label"] == "independent")][["Name","sequence","label","distance"]]
    fullpos = pd.read_csv("output/heterotypic/EtsRunx_v1/seqpositions.csv")[["sequence", "ets_start", "ets_pos", "runx_start", "runx_pos"]]
    dft = dftrain.merge(fullpos, on="sequence").drop_duplicates()
    #dft = dft[dft["distance"] > 4]
    #dft = dft[dft["Name"] == "seq1005_all_clean_seqs_gap0_no20_0.4"]

    # flip pos
    for index, row in dft.iterrows():
        if row['ets_pos'] > row['runx_pos']: # CHANGE THIS DEPENDING ON WHICH CHAMBER
            plus_runx = 1 if row['runx_pos'] - row['runx_start'] == 2 else 2
            dft.at[index,'sequence'] = bio.revcompstr(row["sequence"])
            dft.at[index,'ets_start'] = len(row["sequence"]) - row["ets_start"] - 4
            dft.at[index,'runx_start'] = len(row["sequence"]) - row["runx_start"] - 5
            dft.at[index,'ets_pos'] = dft.at[index,'ets_start'] + 1
            dft.at[index,'runx_pos'] = dft.at[index,'runx_start'] + plus_runx
    print("ets")
    dft["ets_score"], dft["ets_ori"], dft['ets_core'] = pwm_score(dft, pwm_ets, "ets_start", 4, 3)
    print("runx")
    dft["runx_score"], dft["runx_ori"], dft['runx_core'] = pwm_score(dft, pwm_runx, "runx_start", 5, 2)
    #dft['distance'] = dft["runx_pos"] - dft["ets_pos"]

    orimap = {0:"-",1:"+"}
    #dft["orientation"] = dft.apply(lambda x: "%s%s" % (str(x["ets_ori"]), str(x["runx_ori"])),axis=1)
    dft["orientation"] = dft.apply(lambda x: "%s/%s" % (orimap[int(x["ets_ori"])], orimap[int(x["runx_ori"])]),axis=1)
    dft.to_csv("training_pwm.tsv",sep="\t", index=False, float_format='%.3f')

    dft.rename(columns={'runx_score': 'Runx1 strength\n(cooperator TF)', 'ets_score': 'Ets1 strength\n(main TF)'}, inplace=True)
    plot.plot_stacked_categories(dft, "distance", path="distance_bar.png", title="Distance distribution", ratio=True)
    plot.plot_stacked_categories(dft, "orientation", path="ori_bar.png", title="Relative sites orientation\ndistribution", ratio=True)
    plot.plot_box_categories(dft, incols=["Ets1 strength\n(main TF)", "Runx1 strength\n(cooperator TF)"], alternative="smaller")
