import os
os.chdir("../..")
import pandas as pd
import pickle

import chip2probe.modeler.mutation as mut
from chip2probe.probe_generator.probefilter.sitespredict.imads import iMADS
from chip2probe.probe_generator.probefilter.sitespredict.imadsmodel import iMADSModel
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.shapemodel import ShapeModel

if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    df = df[~df['name'].str.contains("dist|weak")]

    imads12_paths = ["input/modeler/imads_model/Ets1_w12_GGAA.model", "input/modeler/imads_model/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19) # 0.2128

    # dm = mut.mutate_dist(df["sequence"].tolist(), imads12, warning=False)
    # dm["label"] = df["label"]
    # dm.to_csv("dm.csv", index=False, header=True)
    dm = pd.read_csv("dm.csv")

    ct = CoopTrain(dm["sequence"].values.tolist(), corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"], imads=imads12)
    ori = ct.get_feature("orientation", {"positive_cores":["GGAA","GGAT"], "relative":True, "one_hot":False})
    dm["orientation"] = [x["ori"] for x in ori]

    feature_dict = {
        "distance":{"type":"numerical"},
         "shape_in": {"seqin":4, "smode":"strength", "direction":"inout"}, # maximum seqin is 4
         "shape_out": {"seqin":-4, "smode":"strength", "direction":"inout"}
    }
    cols = ['dist_numeric', 'ProT_outer_str_pos_2', 'MGW_inner_str_pos_3', 'ProT_outer_wk_pos_2', 'ProT_outer_wk_pos_1', 'Roll_inner_str_pos_2',  'ProT_inner_wk_pos_1', 'HelT_inner_str_pos_2', 'Roll_inner_wk_pos_2', 'MGW_inner_str_pos_2']
    train_df1 = pd.DataFrame(ct.get_feature_all(feature_dict))[cols]
    train1 = train_df1.values.tolist()
    model1 = pickle.load(open("input/modeler/coopmodel/dist_shapeio.sav", "rb"))
    pred1 = model1.predict(train1)
    prob1 = model1.predict_proba(train1)
    dm["shape_pred"] = pred1
    dm["shape_prob"] = [prob1[i][pred1[i]] for i in range(len(pred1))]

    feature_dict = {
        "distance":{"type":"numerical"},
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True},
        "affinity": {"imads":imads12}
    }
    train = ct.get_feature_all(feature_dict, aslist=True)
    model = pickle.load(open("input/modeler/coopmodel/dist_ori_12merimads.sav", "rb"))
    pred = model.predict(train)
    prob = model.predict_proba(train)
    dm["main_pred"] = pred
    dm["main_prob"] = [prob[i][pred[i]] for i in range(len(pred))]
    dm.to_csv("distmut.csv", index=False, header=True)
