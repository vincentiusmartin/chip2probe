import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle

from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.pbmescore import PBMEscore
import chip2probe.modeler.mutation as mut
import chip2probe.training_gen.arranalysis as arr
import chip2probe.training_gen.traingen as tg
from chip2probe.modeler.cooptrain import CoopTrain

def plot_coop(df, lb, suffix, col1="indiv_median", col2="two_median", plot_selected=True):
    dfsub = df[["Name"]].drop_duplicates()
    ax = plt.axes()
    arr.plot_classified_labels(lb[lb["label"] != "anticooperative"], col1=col1, col2=col2, log=True, plotnonsignif=False,
                       xlab="indiv median", ylab="two median", title="Cooperative vs independent binding",
                       labelnames=["cooperative","independent","anticooperative"], axes=ax) #,
    wtavail = lb.merge(dfsub)
    if plot_selected:
        ax.scatter(np.log(wtavail["indiv_median"]), np.log(wtavail["two_median"]), color="cyan", s=1, label="wt_selected")
    ax.legend()
    plt.savefig("in_%s.png" % suffix)
    plt.clf()

def filter_by_delta(df, lb, ncoop=350,nadd=350, uselog=False):
    dfnm = df[["Name"]].drop_duplicates()
    wtavail = lb.merge(dfnm)
    if uselog:
        wtavail["two_median"] = np.log(wtavail["two_median"])
        wtavail["indiv_median"] = np.log(wtavail["indiv_median"])
    wtavail["delta"] = abs(wtavail["two_median"] - wtavail["indiv_median"]) / np.sqrt(2)
    wtcoop = wtavail[wtavail["label"] == "cooperative"].nlargest(ncoop,'delta')[["Name"]]
    wtadd = wtavail[wtavail["label"] == "independent"].nsmallest(nadd,'delta')[["Name"]]
    selected = pd.concat([wtcoop,wtadd])
    return df.merge(selected, on="Name")

def get_rel_ori(ori1,ori2,htth=False):
    if ori1 > 0 and ori2 > 0:
        return "+/+"
    elif ori1 <= 0 and ori2 <= 0:
        return "-/-"
    else:
        if not htth:
            return "+/-" if ori1 > 0 and ori2 <= 0 else "-/+"
        else:
            return "+/-"

def fix_ori(ori):
    return 0 if ori < 0 else ori

def gen_training(seqs, type, kompdict, pwmdict):
    results = []
    if type == "etsets":
        for seq in seqs:
            preds = tg.pred_komp_pwm(seq,kompdict['ets1'],pwmdict['ets1'])
            dictres = {"Sequence":seq,
                       "orientation":get_rel_ori(preds[0]['ori'],preds[1]['ori'],htth=True),
                       "distance":preds[1]["core_mid"] - preds[0]["core_mid"]
                    }
            stridx, wkidx = (0,1) if preds[0]['score'] > preds[1]['score'] else (1,0)
            dictres['site_str_pos'] = preds[stridx]['core_mid']
            dictres['site_str_score'] = preds[stridx]['score']
            dictres['site_str_ori'] = fix_ori(preds[stridx]['ori'])
            dictres['site_wk_pos'] = preds[wkidx]['core_mid']
            dictres['site_wk_score'] = preds[wkidx]['score']
            dictres['site_wk_ori'] = fix_ori(preds[wkidx]['ori'])
            results.append(dictres)
    else: # etsrunx
        for seq in seqs:
            preds = tg.pred_all(seq,kompdict,pwmdict)
            dictres = {"Sequence":seq,
                "orientation":get_rel_ori(preds[0]['ori'],preds[1]['ori'],htth=False), ##
                "distance":preds[1]["core_mid"] - preds[0]["core_mid"]
            }
            eidx, ridx = (0,1) if preds[0]['tf'] == "ets1" else (1,0)
            dictres['ets1_pos'] = preds[eidx]['core_mid']
            dictres['ets1_score'] = preds[eidx]['score']
            dictres['ets1_ori'] = fix_ori(preds[eidx]['ori'])
            dictres['runx1_pos'] = preds[ridx]['core_mid']
            dictres['runx1_score'] = preds[ridx]['score']
            dictres['runx1_ori'] = fix_ori(preds[ridx]['ori'])
            results.append(dictres)
    return results


pd.set_option("display.max_columns",None)
if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"

    # trainingpath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/main_nar/output/Ets1Runx1/training/train_ets1_runx1.tsv"
    # lblpath = "%s/main_nar/output/Ets1Runx1/label_pr/both_ori_plt_ets1_runx1.csv" % basepath
    # train = pd.read_csv(trainingpath,sep="\t")
    # lbled = pd.read_csv(lblpath).rename(columns={"intensity_x":"indiv_median","intensity_y":"two_median"})
    # suffix = "etsrunx"
    # main_modelpath = "%s/main_nar/output/Ets1Runx1/model/ets1_runx1_rfmodel.sav" % basepath
    # shape_modelpath = "%s/main_nar/output/Ets1Runx1/model/rfposmodel.sav" % basepath
    # s1, s2 = "ets1", "runx1"
    # rel_ori = False
    # one_hot_ori = False
    # smode = "positional"

    trainingpath = "%s/main_nar/output/Ets1Ets1_v2/training/train_ets1_ets1.tsv" % basepath
    lblpath = "%s/main_nar/output/Ets1Ets1_v2/label_pr/lbled_o1_selected.csv" % basepath
    train = pd.read_csv(trainingpath,sep="\t")
    lbled = pd.read_csv(lblpath)
    suffix = "etsets"
    main_modelpath = "%s/main_nar/output/Ets1Ets1/model/ets1_ets1_rfmodel.sav" % basepath
    main2_modelpath = "%s/main_nar/output/Ets1Ets1_v2/model/ets1_ets1_rfmodel.sav" % basepath
    shape_modelpath = "%s/main_nar/output/Ets1Ets1/model/rfposmodel.sav" % basepath
    s1, s2 = "site_str", "site_wk"
    rel_ori = True
    one_hot_ori = True
    smode = "relative"

    pwm_ets = PWM("%s/input/sitemodels/pwm/ets1.txt" % basepath, log=True, reverse=False)
    kompas_ets = Kompas("%s/input/sitemodels/kompas/Ets1_kmer_alignment.txt" % basepath, core_start = 11, core_end = 15, core_center = 12)
    escore_ets = PBMEscore("%s/input/sitemodels/escores/Ets1_8mers_11111111.txt" % basepath)

    pwm_runx = PWM("%s/input/sitemodels/pwm/runx1.txt" % basepath, 8, 17, log=True, reverse=True)
    kompas_runx = Kompas("%s/input/sitemodels/kompas/Runx1_kmer_alignment.txt" % basepath, core_start = 12, core_end = 17, core_center = 14)
    escore_runx = PBMEscore("%s/input/sitemodels/escores/Runx1_8mers_11111111.txt" % basepath)

    kompdict = {"ets1":kompas_ets}
    pwmdict = {"ets1":pwm_ets}
    escdict = {"ets1":escore_ets}
    coredict = {"ets1":["GGAA","GGAT"]}
    if suffix == "etsrunx":
        kompdict["runx1"] = kompas_runx
        pwmdict["runx1"] = pwm_runx
        escdict["runx1"] = escore_runx
        coredict["runx1"] = ["GAGGT","GCGGC","GCGGG","GCGGT","GTGGC","GTGGG","GTGGT"]

    plot_coop(train,lbled,"before_%s"%suffix, plot_selected=False)
    train = filter_by_delta(train,lbled,uselog=True).reset_index().head(1000) #ncoop, nadd = 300
    plot_coop(train,lbled,"after_%s"%suffix)

    mutres_all = []
    lbldf = train[["index","label"]].rename(columns={"index":"id","label":"wtlabel"})
    lbldf["wtlabel"] = lbldf["wtlabel"].replace({'cooperative': 1, 'independent': 0})

    wts = mut.initiate_wt(train[["Sequence"]], kompdict, pwmdict)
    mutres_all.extend(wts)

    print("Mutate affinity")
    aff_m = mut.mutate_affinity(train[["Sequence"]], kompdict, pwmdict, coredict)
    mutres_all.extend(aff_m)

    print("Mutate distance")
    aff_d = mut.mutate_distance(train[["Sequence"]], kompdict, pwmdict)
    mutres_all.extend(aff_d)

    print("Mutate orientation")
    aff_o =  mut.mutate_orientation(train[["Sequence"]], kompdict, pwmdict)
    mutres_all.extend(aff_o)

    # TODO: delete wt without mutatants
    custom = pd.DataFrame(mutres_all) \
        .sort_values(["id","comment"]) \
        .merge(lbldf, on="id", how="inner")

    custom_train = gen_training(custom["Sequence"].tolist(), suffix, kompdict, pwmdict)
    ct = CoopTrain(pd.DataFrame(custom_train))

    test = ct.get_feature_all({
        "distance":{"type":"numerical"},
        "affinity": {"colnames": ("%s_score"%s1,"%s_score"%s2)},
        "orientation": {"relative":rel_ori, "one_hot":one_hot_ori, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}}
    })
    mainmodel = pickle.load(open(main_modelpath, "rb"))
    pred = mainmodel.predict(test)
    prob = mainmodel.predict_proba(test)
    custom["main_pred"] = pred
    custom["main_proba"] = [p[1] for p in prob]

    main2model = pickle.load(open(main2_modelpath, "rb"))
    pred2 = main2model.predict(test)
    prob2 = main2model.predict_proba(test)
    custom["main2_pred"] = pred
    custom["main2_proba"] = [p[1] for p in prob]

    shapetest = ct.get_feature_all({
        "distance":{"type":"numerical"},
        "orientation": {"relative":rel_ori, "one_hot":one_hot_ori, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}},
        "shape_in":{"seqin":3, "poscols":['%s_pos'%s1,'%s_pos'%s2], "smode":smode},
        "shape_out":{"seqin":-4, "poscols":['%s_pos'%s1,'%s_pos'%s2], "smode":smode},
        "sequence_in":{"seqin":3, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":smode},
        "sequence_out":{"seqin":-4, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":smode}
    })
    shapemodel = pickle.load(open(shape_modelpath, "rb"))
    shapepred = shapemodel.predict(shapetest)
    shapeprob = shapemodel.predict_proba(shapetest)
    custom["shape_pred"] = shapepred
    custom["shape_proba"] = [p[1] for p in shapeprob]

    custom.to_csv("mutall_%s.csv"%suffix,index=False)
