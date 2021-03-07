import os

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl
from sklearn import ensemble, tree
import subprocess

if __name__ == "__main__":
    # basepath = "output/Ets1Runx1"
    # trainingpath = "%s/training/train_ets1_runx1.tsv" % basepath
    # s1, s2 = "ets1", "runx1"
    # rel_ori = False
    # one_hot_ori = False
    # smode = "positional"

    basepath = "output/Ets1Ets1"
    trainingpath = "%s/training/train_ets1_ets1.tsv" % basepath
    s1, s2 = "site_str", "site_wk"
    rel_ori = True
    one_hot_ori = True
    smode = "relative"

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
                    "orientation": {"relative":rel_ori, "one_hot":one_hot_ori, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}},
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,sequence":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":rel_ori, "one_hot":one_hot_ori, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}},
                    "sequence_in":{"seqin":5, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":smode},
                    "sequence_out":{"seqin":-3, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":smode}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,shape":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":rel_ori, "one_hot":one_hot_ori, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}},
                    "shape_in":{"seqin":5, "poscols":['%s_pos'%s1,'%s_pos'%s2], "smode":smode},
                    "shape_out":{"seqin":-2, "poscols":['%s_pos'%s1,'%s_pos'%s2], "smode":smode}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,sequence,shape":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":rel_ori, "one_hot":one_hot_ori, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}},
                    "shape_in":{"seqin":5, "poscols":['%s_pos'%s1,'%s_pos'%s2], "smode":smode},
                    "shape_out":{"seqin":-2, "poscols":['%s_pos'%s1,'%s_pos'%s2], "smode":smode},
                    "sequence_in":{"seqin":5, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":smode},
                    "sequence_out":{"seqin":-3, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":smode}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("%s_score"%s1,"%s_score"%s2)},
                    "orientation": {"relative":rel_ori, "one_hot":one_hot_ori, "pos_cols": {"%s_pos"%s1:"%s_ori"%s1, "%s_pos"%s2:"%s_ori"%s2}},
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
    }

    pl.plot_model_metrics(best_models, path="%s/model/auc_posfeatures.png"%basepath, cvfold=10, score_type="auc", varyline=True, title="AUC Shape features")
