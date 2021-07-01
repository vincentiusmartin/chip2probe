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
    many_sequences = ["GCCCTGACCCACATCCTGCCGCAGCCTGCTTGTCCT"] #, "CAGCTGGCCGGAACCTGCGTCCCCTTCCCCCGCCGC"]
    df = pd.DataFrame(list(zip(many_sequences, ['seq1','seq2'])), columns=['sequence', 'key'])

    # ========= PWM =========
    pwm_runx = PWM("input/sitemodels/pwm/runx1.txt", 8, 17, log=True, reverse=True)
    pwm_ets = PWM("input/sitemodels/pwm/ets1.txt", log=True, reverse=False)

    # ========= Kompas =========
    kompas_ets = Kompas("input/sitemodels/kompas/Ets1_kmer_alignment.txt", core_start = 11, core_end = 15, core_center = 12)
    kompas_runx = Kompas("input/sitemodels/kompas/Runx1_kmer_alignment.txt", core_start = 12, core_end = 17, core_center = 14)

    kp_ets = KompasPWM(kompas_ets, pwm_ets)
    pred_ets = kp_ets.predict_sequences(many_sequences)
    plot_ets = kp_ets.make_plot_data(pred_ets, color="#b22222")

    kp_runx = KompasPWM(kompas_runx, pwm_runx)
    pred_runx = kp_runx.predict_sequences(many_sequences)
    plot_runx = kp_runx.make_plot_data(pred_runx, color="#0343df")

    # Generate sequence plot
    sp = SitesPlotter()
    sp.plot_seq_combine([plot_ets, plot_runx], filepath="plot.pdf", top_cutoff=10)
