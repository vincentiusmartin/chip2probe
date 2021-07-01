import pandas as pd

import chip2probe.modeler.plotlib as plot

import chip2probe.training_gen.traingen as tg
from sitespredload import imads_ets, imads_runx, kp_ets, kp_runx

def cmp_kp_imads(seq, imads, kompwm, idx=0):
    """
    this assumes only 1 kompas predicted site
    """
    kp_pred = kompwm.predict_sequence(seq)[idx]
    imads_pred = imads.predict_sequence(seq)
    kp_scr, imads_scr = None, None
    # in case of multiple preds
    for p in imads_pred:
        if kp_pred["core_start"] == p["core_start"]:
            kp_scr, imads_scr = kp_pred["score"], p["score"]
            break
    return kp_scr, imads_scr

if __name__ == "__main__":
    fname = "output/mutall_etsets.csv"
    pd.set_option("display.max_columns",None)
    df = pd.read_csv(fname)

    respred = []
    seqlist = list(set(df["Sequence"].tolist()))
    divider = len(seqlist) // 100
    for i in range(0,len(seqlist)):
        if i % divider == 0:
            print("progress %d/%d"%(i,len(seqlist)))
        seq = seqlist[i]
        ets_kp, ets_imads = cmp_kp_imads(seq, imads_ets, kp_ets, 0)
        runx_kp, runx_imads = cmp_kp_imads(seq, imads_ets, kp_ets, 1) ###
        if ets_imads != None and runx_imads != None:
            respred.append({"Sequence":seq, "s1_imads":ets_imads, "s2_imads":runx_imads})
    imads_df = pd.DataFrame(respred)

    df = df.merge(imads_df, on="Sequence").sort_values("id")
    df.to_csv(fname.split(".")[0]+"_with_imads.csv", index=False)
