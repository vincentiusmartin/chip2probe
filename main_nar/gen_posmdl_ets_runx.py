import os

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl
from sklearn import ensemble, tree
import subprocess

if __name__ == "__main__":
    trainingpath = "output/Ets1Runx1/training/train_ets1_runx1.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    ct = CoopTrain(df)
    pd.set_option("display.max_columns",None)

    rf_param_grid = {
        'n_estimators': [500], #[500,750,1000],
        'max_depth': [10], #[5,10,15],
        "min_samples_leaf": [10], #[5,10,15],
        "min_samples_split" : [20], #[5,10,15]
    }

    best_models = {
        "distance,orientation":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets1_pos":"ets1_ori", "runx1_pos":"runx1_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,sequence":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets1_pos":"ets1_ori", "runx1_pos":"runx1_ori"}},
                    "sequence_in":{"seqin":5, "poscols":['ets1_pos','runx1_pos'], "namecol":"Name"},
                    "sequence_out":{"seqin":-3, "poscols":['ets1_pos','runx1_pos'], "namecol":"Name"}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,shape":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets1_pos":"ets1_ori", "runx1_pos":"runx1_ori"}},
                    "shape_in":{"seqin":5, "poscols":['ets1_pos','runx1_pos']},
                    "shape_out":{"seqin":-2, "poscols":['ets1_pos','runx1_pos']}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,sequence,shape":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets1_pos":"ets1_ori", "runx1_pos":"runx1_ori"}},
                    "shape_in":{"seqin":5, "poscols":['ets1_pos','runx1_pos']},
                    "shape_out":{"seqin":-2, "poscols":['ets1_pos','runx1_pos']},
                    "sequence_in":{"seqin":5, "poscols":['ets1_pos','runx1_pos'], "namecol":"Name"},
                    "sequence_out":{"seqin":-3, "poscols":['ets1_pos','runx1_pos'], "namecol":"Name"}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("ets1_score","runx1_score")},
                    "orientation": {"relative":False, "pos_cols": {"ets1_pos":"ets1_ori", "runx1_pos":"runx1_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
    }

    pl.plot_model_metrics(best_models, cvfold=10, score_type="auc", varyline=True, title="Average ROC Curves for Ets1-Runx1\n(using shape and sequence features)")
