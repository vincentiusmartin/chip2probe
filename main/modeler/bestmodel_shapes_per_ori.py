'''
Created on June 9, 2020

@author: Vincentius Martin, Farica Zhuang

Pick the best shape model for different orientation
'''

import pandas as pd
import os
os.chdir("../..")

import chip2probe.modeler.plotlib as pl
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
from chip2probe.modeler.dnashape import DNAShape

if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    # select only genomic (i.e. non-custom) sequences
    df = df[~df['name'].str.contains("dist|weak")]
    cooptr = CoopTrain(df, corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"])
    #
    x_ori = cooptr.get_feature("orientation", {"positive_cores":["GGAA", "GGAT"]})
    df["orientation"] = pd.DataFrame(x_ori)["ori"]

    rf_param_grid = {
        'n_estimators': [500, 1000, 1500],
        'max_depth':[5, 10, 15],
        "min_samples_leaf" : [5, 10, 15],
        "min_samples_split" :[5, 10, 15]
    }

    orientations = ["HH", "TT", "HT/TH"]
    for ori in orientations:
        curdf = df[df["orientation"] == ori]
        curct = CoopTrain(curdf, corelen=4, positive_cores=["GGAA","GGAT"])
        best_models = {
                "distance":
                  BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                              param_grid=rf_param_grid,
                              train_data=curct.get_training_df({
                                      "distance":{"type":"numerical"}
                                  }),
                    ).run_all(),
                "shape":
                  BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                              param_grid=rf_param_grid,
                              train_data=curct.get_training_df({
                                      "shape_in": {"seqin":4, "smode":"positional", "direction":"inout"}, # maximum seqin is 4
                                      "shape_out": {"seqin":-4, "smode":"positional", "direction":"inout"}
                                  }),
                    ).run_all(),
                "dist,shape":
                  BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                              param_grid=rf_param_grid,
                              train_data=curct.get_training_df({
                                      "distance":{"type":"numerical"},
                                      "shape_in": {"seqin":4, "smode":"positional", "direction":"inout"}, # maximum seqin is 4
                                      "shape_out": {"seqin":-4, "smode":"positional", "direction":"inout"}
                                  }),
                    ).run_all(),
                "dist,shape-top10":
                  BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                              param_grid=rf_param_grid,
                              train_data=curct.get_training_df({
                                      "distance":{"type":"numerical"},
                                      "shape_in": {"seqin":4, "smode":"positional", "direction":"inout"}, # maximum seqin is 4
                                      "shape_out": {"seqin":-4, "smode":"positional", "direction":"inout"}
                                  }),
                              topn=10
                    ).run_all(),
        }
        oriname = ori.replace("/","")
        pl.plot_model_metrics(best_models, cvfold=10, score_type="auc", plotname="auc_%s.png" % oriname)
        break
