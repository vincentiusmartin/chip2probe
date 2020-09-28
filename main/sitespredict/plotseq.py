import os
os.chdir("../..")

from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pbmescore import PBMEscore
from chip2probe.sitespredict.kompas import Kompas
import pandas as pd
from chip2probe.sitespredict.sitesplotter import SitesPlotter

if __name__=="__main__":
    single_sequence = "CTGTGTAAATAAAGCCACAGACTTCCGCCTACAGGT"
    many_sequences = ["ACTGGCAGGAAGGGCAGTTTTGGCAGGAAAAGCCAT", "CAGCTGGCCGGAACCTGCGTCCCCTTCCCCCGCCGC"]
    df = pd.DataFrame(list(zip(many_sequences, ['seq1','seq2'])), columns=['sequence', 'key'])

    # ========= PWM =========
    pwm = PWM("input/site_models/pwm/runx1.txt", 7, 18, False)
    pwm_pred = pwm.predict_sequence(single_sequence)
    pwm_pred_list = pwm.predict_sequences(many_sequences)

    # ========= Escore =========
    escore = PBMEscore("input/site_models/escores/Ets1_8mers_11111111.txt")

    escore_pred = escore.predict_sequence(single_sequence)
    escore_pred_list = escore.predict_sequences(many_sequences)
    escore_pred_df = escore.predict_sequences(df, key_colname="key")

    # ========= iMADS =========
    # load imads model, the cores should be loaded in the same order with the model
    imads_paths = ["input/site_models/imads_models/Ets1_w12_GGAA.model", "input/site_models/imads_models/Ets1_w12_GGAT.model"]
    imads_cores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_paths, imads_cores)]
    imads = iMADS(imads_models, 0.19) # 0.2128

    imads_pred = imads.predict_sequence(single_sequence)
    imads_pred_list = imads.predict_sequences(many_sequences)
    imads_pred_df = imads.predict_sequences(df, key_colname="key")

    # ========= Kompas =========
    kompas = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)
    kompas_pred = kompas.predict_sequence(single_sequence)

    # ========= Plot the sequence =========

    # Make the plot objects, make_plot_data accepts prediction result
    pwm_plot = pwm.make_plot_data(pwm_pred_list)
    escore_plot = escore.make_plot_data(escore_pred_list)
    imads_plot = imads.make_plot_data(imads_pred_list)
    kompas_plot = imads.make_plot_data(imads_pred_list)
    # Generate sequence plot
    sp = SitesPlotter()
    sp.plot_seq_combine([imads_plot,escore_plot], filepath="plot.pdf")
