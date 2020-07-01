
from chip2probe.modeler.shapemodel import ShapeModel
import os
import numpy as np
import pandas as pd
from chip2probe.modeler.cooptrain import CoopTrain
from sklearn.metrics import accuracy_score
os.chdir("../..")


if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    df = df[~df['name'].str.contains("dist|weak")]
    df["label"] = df['label'].map({"cooperative":1,"additive":0})
    ct = CoopTrain(df, corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"])
    ori = ct.get_feature("orientation", {"positive_cores":["GGAA", "GGAT"]})
    df["orientation"] = pd.DataFrame(ori)["ori"]
    df = df.loc[~df["orientation"].isnull()]
    print(df.shape)

    top10 = {
        "HH" : {
            "path": "input/modeler/coopmodel/dist_shape_hh.sav",
            "param":['dist_numeric', 'MGW_inner_s2_pos_3', 'Roll_inner_s2_pos_2', 'HelT_inner_s2_pos_2', 'MGW_inner_s1_pos_3', 'ProT_inner_s2_pos_2',  'ProT_inner_s2_pos_3', 'Roll_inner_s1_pos_3', 'ProT_inner_s2_pos_1', 'ProT_inner_s1_pos_1']
            },
        "TT" : {
            "path": "input/modeler/coopmodel/dist_shape_tt.sav",
            "param": ['dist_numeric','Roll_inner_s2_pos_3','HelT_outer_s1_pos_3', 'Roll_inner_s1_pos_0', 'HelT_outer_s2_pos_2', 'ProT_outer_s1_pos_4', 'ProT_outer_s1_pos_2', 'Roll_inner_s2_pos_1', 'MGW_outer_s2_pos_2', 'HelT_inner_s1_pos_2']
            },
        "HT/TH" : {
            "path": "input/modeler/coopmodel/dist_shape_htth.sav",
            "param": ['Roll_inner_s2_pos_2', 'Roll_outer_s2_pos_1', 'ProT_outer_s2_pos_2', 'ProT_inner_s2_pos_1', 'ProT_inner_s2_pos_2', 'Roll_inner_s2_pos_0',  'ProT_inner_s1_pos_1', 'ProT_outer_s2_pos_1', 'ProT_inner_s2_pos_0', 'Roll_inner_s1_pos_0']
            }
    }

    sm = ShapeModel(top10)
    df["pred"], df["proba"] = sm.predict(df)
    print("training accuracy: %.2f", accuracy_score(df["label"],df["pred"]))
