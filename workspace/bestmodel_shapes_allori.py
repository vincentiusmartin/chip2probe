'''
Created on June 9, 2020

@author: Vincentius Martin, Farica Zhuang

Pick the best model
'''

import pandas as pd
import os
os.chdir("..")

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

    rf_param_grid = {
        'n_estimators': [500, 100, 1500],
        'max_depth':[5, 10, 15],
        "min_samples_leaf" : [10, 15, 20],
        "min_samples_split" :[10, 15, 20]
    }

    # TODO: choose per orientation
    best_models = {"dist,ori":
                    BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                          param_grid=rf_param_grid,
                          train_data=cooptr.get_training_df({
                                  "distance":{"type":"numerical"},
                                  "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                              }),
                        topn=10
                ).run_all(),
                "dist,shape_inout":
                    BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                              param_grid=rf_param_grid,
                              train_data=cooptr.get_training_df({
                                      "distance":{"type":"numerical"},
                                      "shape_in": {"seqin":4, "smode":"strength", "direction":"inout"},
                                      "shape_out": {"seqin":-4, "smode":"strength", "direction":"inout"},
                                  }),
                            topn=10
                ).run_all(),
                "dist,shape_ori":
                    BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                              param_grid=rf_param_grid,
                              train_data=cooptr.get_training_df({
                                      "distance":{"type":"numerical"},
                                      "shape_in": {"seqin":4, "smode":"strength", "direction":"orientation", "positive_cores":["GGAA","GGAT"]},
                                      "shape_out": {"seqin":-4, "smode":"strength", "direction":"orientation", "positive_cores":["GGAA","GGAT"]}
                                  }),
                            topn=10
                ).run_all(),
                "dist,ori,shape_inout":
                    BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                              param_grid=rf_param_grid,
                              train_data=cooptr.get_training_df({
                                      "distance":{"type":"numerical"},
                                      "shape_in": {"seqin":4, "smode":"strength", "direction":"inout"},
                                      "shape_out": {"seqin":-4, "smode":"strength", "direction":"inout"},
                                      "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                                  }),
                            topn=10
                ).run_all(),
                "dist,ori,shape_ori":
                    BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                              param_grid=rf_param_grid,
                              train_data=cooptr.get_training_df({
                                      "distance":{"type":"numerical"},
                                      "shape_in": {"seqin":4, "smode":"strength", "direction":"orientation", "positive_cores":["GGAA","GGAT"]},
                                      "shape_out": {"seqin":-4, "smode":"strength", "direction":"orientation", "positive_cores":["GGAA","GGAT"]},
                                      "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                                  }),
                            topn=10
                ).run_all()
    }
    pl.plot_model_metrics(best_models, cvfold=10, score_type="auc", varyline=True)
