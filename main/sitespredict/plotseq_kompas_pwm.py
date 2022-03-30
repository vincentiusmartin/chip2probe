import os
os.chdir("../..")

from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pbmescore import PBMEscore
from chip2probe.sitespredict.kompas import Kompas
import pandas as pd
from chip2probe.sitespredict.sitesplotter import SitesPlotter
from chip2probe.sitespredict.kompaspwm import KompasPWM

import matplotlib.pyplot as plt
import matplotlib.patches as patches

if __name__=="__main__":
    # many_sequences = ["AATCAGAGGAAACAGATAGACTGAGGAAGTGAAGAA","AATCAGAGGAAACAGATAGACTGAGGAGGTGAAGAA","AATCAGAGGACACAGATAGACTGAGGAAGTGAAGAA","AATCAGAGGACACAGATAGACTGAGGAGGTGAAGAA","GATCCCAAACAGGATATCTGTGGTAAGCA"] #, "CAGCTGGCCGGAACCTGCGTCCCCTTCCCCCGCCGC"]
    # df = pd.DataFrame(list(zip(many_sequences, ['wt','m1','m2','m3'])), columns=['sequence', 'key'])

    df = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/main_nar/input/probefiles/Ets1_Runx1_pr_clean.csv")
    df = df[(df["ori"] == "o1") & (df["rep"] == "r1")].rename(columns={"Sequence": "sequence", "Name": "key"})[["sequence","key"]].head(100)

    # ========= PWM =========
    pwm_runx = PWM("input/sitemodels/pwm/runx1.txt", 8, 17, log=True, reverse=True)
    pwm_ets = PWM("input/sitemodels/pwm/ets1.txt", log=True, reverse=False)

    # ========= Kompas =========
    kompas_ets = Kompas("input/sitemodels/kompas/Ets1_kmer_alignment.txt", core_start = 11, core_end = 15, core_center = 12)
    kompas_runx = Kompas("input/sitemodels/kompas/Runx1_kmer_alignment.txt", core_start = 12, core_end = 17, core_center = 14)

    kp_ets = KompasPWM(kompas_ets, pwm_ets)
    pred_ets = kp_ets.predict_sequences(df)
    plot_ets = kp_ets.make_plot_data(pred_ets, color="#b22222")

    kp_runx = KompasPWM(kompas_runx, pwm_runx)
    pred_runx = kp_runx.predict_sequences(df)
    plot_runx = kp_runx.make_plot_data(pred_runx, color="#0343df")

    # Generate sequence plot
    sp = SitesPlotter()
    sp.plot_seq_combine([plot_ets, plot_runx], filepath="plot.pdf", top_cutoff=8)
