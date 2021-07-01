import pandas as pd
import pickle

import chip2probe.training_gen.traingen as tg
from make_custom import gen_training
import chip2probe.modeler.mutation as mut

from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.kompas import Kompas

from chip2probe.modeler.cooptrain import CoopTrain

basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"

# path = "input/litseqs_etsrunx.csv"
# suffix = "etsrunx"
# s1, s2 = "ets1", "runx1"
# rel_ori = False
# one_hot_ori = False
# smode = "positional"
# main_modelpath = "%s/main_nar/output/Ets1Runx1/model/ets1_runx1_rfmodel.sav" % basepath
# shape_modelpath = "%s/main_nar/output/Ets1Runx1/model/rfposmodel.sav" % basepath

path = "input/litseqs_etsets.csv"
suffix = "etsets"
s1, s2 = "site_str", "site_wk"
rel_ori = True
one_hot_ori = True
smode = "relative"
main_modelpath = "%s/main_nar/output/Ets1Ets1_v2/model/ets1_ets1_rfmodel.sav" % basepath
shape_modelpath = "%s/main_nar/output/Ets1Ets1/model/rfposmodel.sav" % basepath

pwm_runx = PWM("%s/input/sitemodels/pwm/runx1.txt"%basepath, 8, 17, log=True, reverse=True)
pwm_ets = PWM("%s/input/sitemodels/pwm/ets1.txt"%basepath, log=True, reverse=False)
kompas_ets = Kompas("%s/input/sitemodels/kompas/Ets1_kmer_alignment.txt"%basepath, core_start = 11, core_end = 15, core_center = 12)
kompas_runx = Kompas("%s/input/sitemodels/kompas/Runx1_kmer_alignment.txt"%basepath, core_start = 12, core_end = 17, core_center = 14)

kompdict = {"ets1":kompas_ets}
pwmdict = {"ets1":pwm_ets}
if suffix == "etsrunx":
    kompdict["runx1"] = kompas_runx
    pwmdict["runx1"] = pwm_runx

df = pd.read_csv(path)
seqs = df["Sequence"].tolist()

res = []
for idx,row in df.iterrows():
    pred = tg.pred_all(row["Sequence"], kompdict, pwmdict)
    res.append(mut.makeseqdict(row["Sequence"], -1, pred, type="literature", comment=row["id"]))
resdf = pd.DataFrame(res)

train = gen_training(seqs, suffix, kompdict, pwmdict)
ct = CoopTrain(pd.DataFrame(train))

test = ct.get_feature_all({
    "distance":{"type":"numerical"},
    "affinity": {"colnames": ("%s_score"%s1,"%s_score"%s2)},
    "orientation": {"relative":rel_ori, "one_hot":one_hot_ori, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}}
})
mainmodel = pickle.load(open(main_modelpath, "rb"))
pred = mainmodel.predict(test)
prob = mainmodel.predict_proba(test)
resdf["main_pred"] = pred
resdf["main_proba"] = [p[1] for p in prob]

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
resdf["shape_pred"] = shapepred
resdf["shape_proba"] = [p[1] for p in shapeprob]

resdf["select"] = ""

resdf.to_csv("sitedata_%s.csv"%suffix,index=False)
