import pandas as pd

from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.kompaspwm import KompasPWM

import matplotlib.pyplot as plt
import numpy as np

pd.set_option("display.max_columns",None)

def cmp_kp_imads(seq, imads, kompwm):
    """
    this assumes only 1 kompas predicted site
    """
    kp_pred = kompwm.predict_sequence(seq)[0]
    imads_pred = imads.predict_sequence(seq)
    kp_scr, imads_scr = None, None
    # in case of multiple preds
    for p in imads_pred:
        if kp_pred["core_start"] == p["core_start"]:
            kp_scr, imads_scr = kp_pred["score"], p["score"]
            break
    return kp_scr, imads_scr

def checkdelta(imads, kompwm, df):
    grps = df.groupby("id")
    imads_deltas, kp_deltas = [], []
    for idx, group in grps:
        if group[group["comment"] == "wt"].empty:
            print("empty",idx)
            continue
        wtseq = group[group["comment"] == "wt"]["Sequence"].tolist()[0]
        wtkp, wtimads = cmp_kp_imads(wtseq, imads, kompwm)
        if wtkp == None or wtimads == None:
            continue
        otherseqs = group[group["comment"] != "wt"]["Sequence"].tolist()
        for seq in otherseqs:
            scrkp, scrimads = cmp_kp_imads(seq, imads, kompwm)
            if scrkp != None or scrimads != None:
                imads_deltas.append(scrimads-wtimads)
                kp_deltas.append(scrkp-wtkp)
    return imads_deltas, kp_deltas

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
imads_runx = iMADS(imads_runx_models, 0.25) # 0.2128

# ========= PWM =========
pwm_runx = PWM("%s/input/sitemodels/pwm/runx1.txt"%basesp, 8, 17, log=True, reverse=True)
pwm_ets = PWM("%s/input/sitemodels/pwm/ets1.txt"% basesp, log=True, reverse=False)

# ========= KOMPAS =========
kompas_ets = Kompas("%s/input/sitemodels/kompas/Ets1_kmer_alignment.txt"% basesp, core_start = 11, core_end = 15, core_center = 12)
kompas_runx = Kompas("%s/input/sitemodels/kompas/Runx1_kmer_alignment.txt"% basesp, core_start = 12, core_end = 17, core_center = 14)

# ========= KompasPWM =========
kp_ets = KompasPWM(kompas_ets, pwm_ets)
kp_runx = KompasPWM(kompas_runx, pwm_runx)

df = pd.read_csv("output/selected_er_fin.csv")
df = df[(df["comment"] == "wt") | (df['comment'].str.contains("s1") & (df["muttype"] == "strength"))]
imads_deltas, kp_deltas = checkdelta(imads_ets, kp_ets, df)

cordf = pd.DataFrame(list(zip(imads_deltas,kp_deltas)), columns=["imads_delta", "pwm_delta"])
x, y = cordf["imads_delta"].tolist(), cordf["pwm_delta"].tolist()
r_squared = np.corrcoef(x,y)[0,1]**2
print(r_squared)
plt.plot(x,y,color="blue")
cordf.plot.scatter(x='imads_delta', y='pwm_delta', c='Blue', s=0.8)

slope, intercept = np.polyfit(x, y, 1)
abline_values = [slope * i + intercept for i in x]
plt.plot([min(x),max(x)], [min(abline_values),max(abline_values)], color='red', linestyle='dashed')
# plt.text(-7.5, 0, 'R²=%.2f'%r_squared, size=12, color='red')

plt.title("Delta ets1")
plt.savefig('scatter_ets1.png')
