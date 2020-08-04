from chip2probe.modeler.shapemodel import ShapeModel
from chip2probe.probe_generator.probefilter.sitespredict.imads import iMADS
from chip2probe.probe_generator.probefilter.sitespredict.imadsmodel import iMADSModel
import os
import numpy as np
import pandas as pd
from chip2probe.modeler.cooptrain import CoopTrain
from sklearn.metrics import accuracy_score
import chip2probe.util.coopgeneral as cg
os.chdir("../..")

if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    poscores = ["GGAA", "GGAT"]

    imads12_paths = ["input/modeler/imads_model/Ets1_w12_GGAA.model", "input/modeler/imads_model/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19)

    df = pd.read_csv(trainingpath, sep="\t")
    df = df[~df['name'].str.contains("dist|weak")]
    df["label"] = df['label'].map({"cooperative":1,"additive":0})
    df["orientation"] = df.apply(lambda x: cg.get_relative_orientation(x["sequence"],imads12), axis = 1)

    top10 = {
        "HH" : {
            "path": "input/modeler/coopmodel/dist_shape_hh.sav",
            "param":['Roll_inner_s2_pos_2', 'Roll_outer_s2_pos_1', 'ProT_outer_s2_pos_2', 'ProT_inner_s2_pos_1', 'ProT_inner_s2_pos_2', 'Roll_inner_s2_pos_0',  'ProT_inner_s1_pos_1', 'ProT_outer_s2_pos_1', 'ProT_inner_s2_pos_0', 'Roll_inner_s1_pos_0']
            },
        "TT" : {
            "path": "input/modeler/coopmodel/dist_shape_tt.sav",
            "param":['dist_numeric','Roll_inner_s2_pos_3','HelT_outer_s1_pos_3', 'Roll_inner_s1_pos_0', 'HelT_outer_s2_pos_2', 'ProT_outer_s1_pos_4', 'ProT_outer_s1_pos_2', 'Roll_inner_s2_pos_1', 'MGW_outer_s2_pos_2', 'HelT_inner_s1_pos_2']
            },
        "HT/TH" : {
            "path": "input/modeler/coopmodel/dist_shape_htth.sav",
            "param":  ['dist_numeric', 'MGW_inner_s2_pos_3', 'Roll_inner_s2_pos_2', 'HelT_inner_s2_pos_2', 'MGW_inner_s1_pos_3', 'ProT_inner_s2_pos_2',  'ProT_inner_s2_pos_3', 'Roll_inner_s1_pos_3', 'ProT_inner_s2_pos_1', 'ProT_inner_s1_pos_1']
            }
    }

    sm = ShapeModel(top10)
    df["pred"], df["proba"] = sm.predict(df)
    print("training accuracy: %.4f" % accuracy_score(df["label"],df["pred"]))
