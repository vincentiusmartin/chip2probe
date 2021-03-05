import pandas as pd
import os
import pickle
import sklearn.metrics as met
from chip2probe.modeler.cooptrain import CoopTrain
import chip2probe.training_gen.arranalysis as arr
from scipy.stats import pearsonr
import chip2probe.modeler.plotlib as pl
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt
os.chdir("../../../..")

# # get the selected wt
# wtseqs = selected[selected["comment"] == "wt"][["sequence"]].drop_duplicates()
# selwtids = seqlbled.merge(wtseqs, on="sequence")[["id"]]

def make_label(po1, po2, coop_pcut):
    if po1 < coop_pcut and po2 < coop_pcut:
        return 1
    else:
        return 0


if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    selected = pd.read_csv("output/array_design_files/Coop2Ets_validation/custom_probes_selected_repred.csv")
    sl = pd.read_csv("output/homotypic/custom_sequences/20nM/05_nofdr_either_ori/seqlbled.csv")

    sl = sl[(sl["label_o1"] != "anticoop") & (sl["label_o2"] != "anticoop")]
    coop_pcut = 0.05
    add_pcut = 0.4
    sl = sl[((sl["p_o1"] < coop_pcut) &  (sl["p_o2"] < coop_pcut)) |
            ((sl["p_o1"] > add_pcut) & (sl["p_o2"] > add_pcut) & (sl["label_o1"] == "additive") & (sl["label_o2"] == "additive"))]
    sl['label'] = sl.apply(lambda x: make_label(x["p_o1"],x["p_o2"],coop_pcut),axis=1)

    # get the sequence ids from the array ids
    orig_df = pd.read_csv("output/array_design_files/Coop2Ets_validation/arrayseqs.csv")
    orig_df = orig_df.loc[orig_df["id"].str.contains("_wt")][["id"]]
    allorig = sl.merge(orig_df, on=["id"])

    alldf = allorig.merge(selected, on="sequence")
    alldf = alldf.drop_duplicates(subset=['sequence'])
    alldf = alldf.sort_values(["muttype","seqid","distance"], ascending=[True,True,False])
    alldf.to_csv("05coop_6add.csv",index=False)
    predunique = alldf[["sequence","id","label","main_pred"]].drop_duplicates()
    alladd = predunique[predunique["label"] == 0]
    matchadd = alladd[alladd["main_pred"] == 0]
    percentadd = float(matchadd.shape[0]) / alladd.shape[0] * 100
    allcoop = predunique[predunique["label"] == 1]
    matchcoop = allcoop[allcoop["main_pred"] == 1]
    percentcoop = float(matchcoop.shape[0]) / allcoop.shape[0] * 100

    print("Matching add %d/%d (%.2f%%), coop %d/%d (%.2f%%)" % (matchadd.shape[0], alladd.shape[0], percentadd, matchcoop.shape[0], allcoop.shape[0], percentcoop), "\n")

    g = alldf.groupby("muttype")

    totalwt = alldf[(alldf["comment"] == "wt")]
    wtadd = alldf[(alldf["comment"] == "wt") & (alldf["label"] == 0)].shape[0]
    wtadd_pred = alldf[(alldf["comment"] == "wt") & (alldf["main_pred"] == 0)].shape[0]
    wtcoop = alldf[(alldf["comment"] == "wt") & (alldf["label"] == 1)].shape[0]
    wtcoop_pred = alldf[(alldf["comment"] == "wt") & (alldf["main_pred"] == 1)].shape[0]
    print("add",wtadd,wtadd_pred,"coop",wtcoop,wtcoop_pred)

    fpr, tpr, _ = metrics.roc_curve(np.array(alldf["label"]), np.array(alldf["main_prob"]))
    auc = metrics.auc(fpr, tpr)
    plt.plot(fpr,tpr, label='AUC all = %.2f' % auc)

    for mtype in ["distance","affinity","orientation"]:
        print("======================================")
        print(mtype)
        curdf = alldf[alldf["muttype"] == mtype]
        curdf_add = curdf[curdf["label"] == 0]
        pred_add = curdf_add[curdf_add["main_pred"] == 0].shape[0]
        addfrac = float(pred_add)/curdf_add.shape[0]*100
        curdf_coop = curdf[curdf["label"] == 1]
        pred_coop = curdf_coop[curdf_coop["main_pred"] == 1].shape[0]
        coopfrac = float(pred_coop)/curdf_coop.shape[0]*100
        cm = met.confusion_matrix(curdf["label"], curdf["main_pred"],labels=[1,0])
        print("Total correct predictions: add %d/%d (%.2f%%), coop %d/%d (%.2f%%)" % (pred_add, curdf_add.shape[0], addfrac, pred_coop, curdf_coop.shape[0], coopfrac))
        print("Confusion matrix:")
        print(cm)
        minp = curdf.apply(lambda x: x["p_o1"] if x["p_o1"]  < x["p_o2"] else x["p_o2"],axis=1)
        corr, _ = pearsonr(minp, curdf["main_prob"])
        print("R(model_probability,p_value): %.2f" %   corr,"\n")
        fprtype, tprtype, _ = metrics.roc_curve(np.array(curdf["label"]), np.array(curdf["main_prob"]))
        auctype = metrics.auc(fprtype, tprtype)
        plt.plot(fprtype,tprtype, label='AUC %s = %.2f' % (mtype,auctype))
        # if mtype == "distance":
        #     g = curdf.groupby(["distance","label"])["label"].count()
        #     pl.plot_stacked_categories(curdf, "distance", path="cust_stackedbar_dist.png", ratio=True)
        # elif mtype == "affinity":
        #     g = curdf.groupby("seqid")
        #     ghead = g.head(1)
        #     ghead = ghead[ghead["label"] == 0][["seqid","id"]]
        #     gtail = g.tail(1)
        #     gtail = gtail[gtail["label"] == 1][["seqid","id"]]
        #     ght = ghead.merge(gtail, on="seqid")
        #     print(ght)
        #     # if needed, print
        # elif mtype == "orientation":
        #     for ori in ["HH","HT/TH","TT"]:
        #         oridf = curdf[(curdf["orientation"] == ori)]
        #         orimatch = oridf[oridf["label"] == oridf["main_pred"]].shape[0]
        #         print("Orientation %s, match: (%d/%d)" % (ori,orimatch,oridf.shape[0]))
        #     #pl.plot_stacked_categories(curdf, "orientation", path="cust_stackedbar_ori.png", ratio=True)

    plt.legend(loc='lower right')
    plt.title("ROC curve for validation array")
    plt.savefig("auc.png")
    plt.clf()
