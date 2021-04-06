import pandas as pd
import pickle
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt

from traingen_ets_ets import gen_training
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.kompas import Kompas

def fill_arr_data(df):
    # we need to map this back to the original array since the custom array's sequences had some modifications
    arrseqs = pd.read_csv("input/ets1ets1validation/arrseqs.csv")
    arrseqs = arrseqs.loc[arrseqs["id"].str.contains("_wt")][["id"]].rename(columns={"id":"Name"})
    arrseqs['Name'] = arrseqs['Name'].str.split('_').str[0]

    arrdata = pd.read_csv("input/ets1ets1validation/arrdata.csv")

    arrdata_sub = arrdata[["Sequence","comment"]]
    arrdata_sub = arrdata_sub[arrdata_sub["comment"] == "wt"].drop_duplicates()
    print(arrdata_sub.count())

    dfinfo = df.merge(arrseqs, on="Name").merge(arrdata, on="Sequence")
    return dfinfo[["Name","Sequence","label","wtlabel","comment"]]

def make_label(po1, po2, coop_pcut):
    if po1 < coop_pcut and po2 < coop_pcut:
        return 1
    else:
        return 0

def print_matching(df, comparecol):
    allindep = df[df["label"] == 0]
    matchindep = allindep[allindep[comparecol] == 0]
    percentindep = float(matchindep.shape[0]) / allindep.shape[0] * 100

    allcoop = df[df["label"] == 1]
    matchcoop = allcoop[allcoop[comparecol] == 1]
    percentcoop = float(matchcoop.shape[0]) / allcoop.shape[0] * 100
    print("Matching on %s:\n indep %d/%d (%.2f%%)\n coop %d/%d (%.2f%%)" % (comparecol, matchindep.shape[0], allindep.shape[0], percentindep, matchcoop.shape[0], allcoop.shape[0], percentcoop), "\n")

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    df = pd.read_csv("output/Ets1Ets1_validation/label_pr/ets_ets_seqlabeled.csv").drop_duplicates()
    df = df[(df["label_o1"] != "anticoop") & (df["label_o2"] != "anticoop")] # filter out anticoop
    df["Name"] = df["Name"].astype(str)

    coop_pcut = 0.1
    indep_pcut = 0.8
    # filter based on: either < pcoop both orientations or independent in both orientations with p val larger than cutoff
    df = df[((df["p_o1"] <= coop_pcut) &  (df["p_o2"] <= coop_pcut)) |
            ((df["p_o1"] > indep_pcut) & (df["p_o2"] > indep_pcut) & (df["label_o1"] == "independent") & (df["label_o2"] == "independent"))]
    df['label'] = df.apply(lambda x: make_label(x["p_o1"],x["p_o2"],coop_pcut),axis=1)

    print("Label count:\n", df["label"].value_counts())

    df = fill_arr_data(df)
    print("------------")
    print_matching(df, "wtlabel")

    # genrate training data
    pwm_ets = PWM("input/sitemodels/ets1.txt", log=True)
    kompas_ets = Kompas("input/sitemodels/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)
    train = gen_training(df, pwm_ets, kompas_ets)
    train.to_csv("output/Ets1Ets1_validation/training/train_ets_val.csv",index=False)

    """
    ct = CoopTrain(train)
    # Predict the sequence using dist, orientation, strength
    feature_dict = {
        "distance":{"type":"numerical"},
        "affinity": {"colnames": ("site_str_score","site_wk_score")},
        "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
    }
    model = pickle.load(open("output/Ets1Ets1/model/ets1_ets1_rfmodel.sav", "rb"))
    trainlist = pd.DataFrame(ct.get_feature_all(feature_dict)).values.tolist()

    pred = model.predict(trainlist)
    prob = model.predict_proba(trainlist)
    train["main_pred"] = pred
    train["main_prob"] = [prob[i][1] for i in range(len(prob))]

    # Predict the sequence using shape
    s1, s2 = "site_str", "site_wk"
    feature_dict = {
        "distance":{"type":"numerical"},
        "orientation": {"relative":True, "one_hot":True, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}},
        "shape_in":{"seqin":5, "poscols":['%s_pos'%s1,'%s_pos'%s2], "smode":"relative"},
        "shape_out":{"seqin":-2, "poscols":['%s_pos'%s1,'%s_pos'%s2], "smode":"relative"}
    }
    shapemodel = pickle.load(open("output/Ets1Ets1/model/ets1_ets1_rfshapemodel.sav", "rb"))
    trainlist = pd.DataFrame(ct.get_feature_all(feature_dict)).values.tolist()

    shapepred = shapemodel.predict(trainlist)
    shapeprob = shapemodel.predict_proba(trainlist)
    train["shape_pred"] = shapepred
    train["shape_prob"] = [shapeprob[i][1] for i in range(len(shapeprob))]

    train.to_csv("output/Ets1Ets1_validation/training/train_wpred.csv",index=False)
    """

    train = pd.read_csv("output/Ets1Ets1_validation/training/train_wpred.csv")
    for p in [("main_prob", "distance,orientation,strength"),("shape_prob", "distance,orientation,shape")]:
        probcol, legend = p
        fpr, tpr, _ = metrics.roc_curve(np.array(train["label"]), np.array(train[probcol]))
        auc = metrics.auc(fpr, tpr)
        plt.plot(fpr,tpr, label='%s = %.2f' % (legend,auc))
    leg = plt.legend(loc="lower right",title="AUC for each combination:")
    leg._legend_box.align = "left"
    plt.savefig("auc.png")
