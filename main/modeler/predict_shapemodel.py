
from chip2probe.modeler.shapemodel import ShapeModel
import os
import pandas as pd
os.chdir("../..")


if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    df = df[~df['name'].str.contains("dist|weak")]

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
    sm.predict(df)
