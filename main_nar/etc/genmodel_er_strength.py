import os
os.chdir("..")

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl
import pickle
from sklearn import ensemble, tree
import subprocess


if __name__ == "__main__":
    basepath = "output/Ets1Runx1"
    trainingpath = "output/Ets1Runx1/training/train_ets1_runx1.tsv"

    df = pd.read_csv(trainingpath, sep="\t")
    ct = CoopTrain(df)
    pd.set_option("display.max_columns",None)

    rf_param_grid = {
        'n_estimators': [1000], #[500,750,1000],
        'max_depth': [10], #[5,10,15],
        "min_samples_leaf": [5], #[5,10,15],
        "min_samples_split" : [5] #[5,10,15]
    }

    best_models = {
        "Runx1 (cooperator TF) strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "affinity": {"colnames": ["runx1_score"]}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "Ets1 (main TF) strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "affinity": {"colnames": ["ets1_score"]}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,Runx1 strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ["runx1_score"]},
                    "orientation": {"relative":False, "pos_cols": {"ets1_pos":"ets1_ori", "runx1_pos":"runx1_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,Ets1 strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ["ets1_score"]},
                    "orientation": {"relative":False, "pos_cols": {"ets1_pos":"ets1_ori", "runx1_pos":"runx1_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "strength,distance,orientation":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("ets1_score","runx1_score")},
                    "orientation": {"relative":False, "pos_cols": {"ets1_pos":"ets1_ori", "runx1_pos":"runx1_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
    }

    pl.plot_model_metrics(best_models, path="%s/model/auc_separate_str.png" % basepath, cvfold=10, score_type="auc", varyline=True, title="Average ROC Curves for Ets1-Runx1", interp=True)

    # tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
    #         feature_names = train.columns,
    #         class_names = ['additive','cooperative'],
    #         rounded = True, proportion = False,
    #         precision = 2, filled = True)
    # subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])
