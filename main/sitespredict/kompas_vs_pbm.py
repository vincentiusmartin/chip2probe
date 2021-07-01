import os
os.chdir("../..")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.kompaspwm import KompasPWM
from chip2probe.sitespredict.pbmescore import PBMEscore

from chip2probe.sitespredict.sitesplotter import SitesPlotter

def plotcor(imads, kompwm, df, tfname):
    lcor1, lcor2 = [], []
    for i, row in df.iterrows():
        cur_imads_scr, cur_kp_scr = None, None
        for imadspred in imads[i]:
            if imadspred["core_start"] == row["%s_start" % tfname]:
                cur_imads_scr = imadspred["score"]
                break
        for kp in kompwm[i]:
            if kp["core_start"] == row["%s_start" % tfname]:
                cur_kp_scr = kp["score"]
                break
        if cur_imads_scr != None and cur_kp_scr != None:
            lcor1.append(cur_imads_scr)
            lcor2.append(cur_kp_scr)
    cordf = pd.DataFrame(list(zip(lcor1,lcor2)), columns=["imads_score", "pwm_score"])
    x, y = cordf["imads_score"].tolist(), cordf["pwm_score"].tolist()
    print("Rsq", cordf["imads_score"].corr(cordf["pwm_score"])**2)
    plt.plot(x,y,color="blue")
    cordf.plot.scatter(x='imads_score', y='pwm_score', c='Blue', s=0.8)

    slope, intercept = np.polyfit(x, y, 1)
    abline_values = [slope * i + intercept for i in x]
    plt.plot([min(x),max(x)], [min(abline_values),max(abline_values)], color='red', linestyle='dashed')

    plt.savefig('scatter_%s.png'%tfname)

if __name__ == "__main__":
    seqdf = pd.read_csv("main_nar/output/Ets1Runx1/training/train_ets1_runx1.tsv", sep="\t").drop_duplicates("Sequence").reset_index()
    seqlist = seqdf["Sequence"].tolist()
    # seqlist = ["GGGGCAGCTTCCGGCAGGAGAGACCACATTCACGGC", "GGGGCAGCTTCTGGCAGGAGAGACCACATTCACGGC", "GGGGCAGCTTCCGGCAGGAGAGACCAGATTCACGGC"]

    # ========= iMADS =========
    imads_ets_paths = ["input/sitemodels/imads_models/Ets1_w12_GGAA.model", "input/sitemodels/imads_models/Ets1_w12_GGAT.model"]
    imads_ets_cores = ["GGAA", "GGAT"]
    imads_ets_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_ets_paths, imads_ets_cores)]
    imads_ets = iMADS(imads_ets_models, 0.19) # 0.2128

    imads_runx_paths = ["input/sitemodels/imads_models/Runx1_w20_GAGGT.model",
                        "input/sitemodels/imads_models/Runx1_w20_GCGGC.model",
                        "input/sitemodels/imads_models/Runx1_w20_GCGGG.model",
                        "input/sitemodels/imads_models/Runx1_w20_GCGGT.model",
                        "input/sitemodels/imads_models/Runx1_w20_GTGGC.model",
                        "input/sitemodels/imads_models/Runx1_w20_GTGGG.model",
                        "input/sitemodels/imads_models/Runx1_w20_GTGGT.model"]
    imads_runx_cores = ["GAGGT", "GCGGC", "GCGGG", "GCGGT", "GTGGC", "GTGGG", "GTGGT"]
    imads_runx_models = [iMADSModel(path, core, 20, [1, 2, 3]) for path, core in zip(imads_runx_paths, imads_runx_cores)]
    imads_runx = iMADS(imads_runx_models, 0.25) # 0.2128

    pred_imads_ets = imads_ets.predict_sequences(seqlist)
    plot_imads_ets = imads_ets.make_plot_data(pred_imads_ets,color="#FFA07A")
    pred_imads_runx = imads_runx.predict_sequences(seqlist)
    plot_imads_runx = imads_runx.make_plot_data(pred_imads_runx, color="#75bbfd")

    # ========= PWM =========
    pwm_runx = PWM("input/sitemodels/pwm/runx1_jaspar.txt", log=True, reverse=False) # 8, 17,
    pwm_ets = PWM("input/sitemodels/pwm/ets1.txt", log=True, reverse=False)
    pred_pwm_ets = pwm_ets.predict_sequences(seqlist)
    pred_pwm_runx = pwm_runx.predict_sequences(seqlist)
    plot_pwm_ets = pwm_ets.make_plot_data(pred_pwm_ets, color="#FFA07A")
    plot_pwm_runx = pwm_runx.make_plot_data(pred_pwm_runx, color="#75bbfd")

    # ========= E-score =========
    escore_ets = PBMEscore("input/sitemodels/escores/Ets1_8mers_11111111.txt")
    escore_runx = PBMEscore("input/sitemodels/escores/Runx1_8mers_11111111.txt")

    pred_escore_ets = escore_ets.predict_sequences(seqlist)
    plot_escore_ets = escore_ets.make_plot_data(pred_escore_ets,color="red")
    pred_escore_runx = escore_runx.predict_sequences(seqlist)
    plot_escore_runx = escore_runx.make_plot_data(pred_escore_runx,color="blue")


    # ========= KOMPAS =========
    kompas_ets = Kompas("input/sitemodels/kompas/Ets1_kmer_alignment.txt", core_start = 11, core_end = 15, core_center = 12)
    kompas_runx = Kompas("input/sitemodels/kompas/Runx1_kmer_alignment.txt", core_start = 12, core_end = 17, core_center = 14)

    pred_kompas_ets = kompas_ets.predict_sequences(seqlist)
    pred_kompas_runx = kompas_runx.predict_sequences(seqlist)
    plot_kompas_ets = kompas_ets.make_plot_data(pred_kompas_ets, color="#FFA07A")
    plot_kompas_runx = kompas_runx.make_plot_data(pred_kompas_runx, color="#75bbfd")

    # ========= KompasPWM =========
    kp_ets = KompasPWM(kompas_ets, pwm_ets)
    kp_runx = KompasPWM(kompas_runx, pwm_runx)

    pred_kp_ets = kp_ets.predict_sequences(seqlist)
    plot_kp_ets = kp_ets.make_plot_data(pred_kp_ets, color="#FFA07A")
    pred_kp_runx = kp_runx.predict_sequences(seqlist)
    plot_kp_runx = kp_runx.make_plot_data(pred_kp_runx, color="#75bbfd")

    l1 = [pred_imads_runx[k].predictions for k in pred_imads_runx]
    l2 = [pred_kp_runx[k].predictions for k in plot_kp_runx]
    print(len(l1),len(l2))
    plotcor(l1,l2,seqdf,"runx1")

    l1 = [pred_imads_ets[k].predictions for k in pred_imads_ets]
    l2 = [pred_kp_ets[k].predictions for k in plot_kp_ets]
    print(len(l1),len(l2))
    plotcor(l1,l2,seqdf,"ets1")

    # ========= Generate sequence plot =========
    sp = SitesPlotter()
    sp.plot_seq_combine([plot_pwm_ets,plot_pwm_runx], filepath="plot_pwm.pdf",top_cutoff=12)
    sp.plot_seq_combine([plot_kompas_ets,plot_kompas_runx], filepath="plot_kompas.pdf")
    sp.plot_seq_combine([plot_imads_ets,plot_imads_runx,plot_escore_ets,plot_escore_runx], filepath="plot_imads.pdf")
    sp.plot_seq_combine([plot_kp_ets,plot_kp_runx], filepath="plot_kp.pdf", top_cutoff=12)
