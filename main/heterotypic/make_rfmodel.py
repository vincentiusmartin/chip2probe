import os
os.chdir("../..")

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl
from sklearn import ensemble, tree
import subprocess

if __name__ == "__main__":
    trainingpath = "output/heterotypic/EtsRunx_v1/ch1_ch2/training_pwm.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    ct = CoopTrain(df)
    pd.set_option("display.max_columns",None)

    rf_param_grid = {
        'n_estimators': [500],
        'max_depth':[5],
        "min_samples_leaf" : [10],
        "min_samples_split" :[10]
    }

    best_models = {
        "distance":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"}
                })
            ).run_all(),
        "orientation":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}
                })
            ).run_all(),
        "strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "affinity": {"colnames": ("ets_score","runx_score")}
                })
            ).run_all(),
        "distace,orientation":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}
                })
            ).run_all(),
        "distance,strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("ets_score","runx_score")}
                })
            ).run_all(),
        "distance,orientation,strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("ets_score","runx_score")},
                    "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}
                })
            ).run_all(),
    }
    pl.plot_model_metrics(best_models, cvfold=10, score_type="auc", varyline=True)

    feature_dict = {
        "distance":{"type":"numerical"},
        "affinity": {"colnames": ("ets_score","runx_score")},
        "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}
    }
    train = ct.get_feature_all(feature_dict)
    label = ct.get_numeric_label({'cooperative': 1, 'additive': 0})
    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=5, min_samples_leaf=10, min_samples_split=10)
    m = rf.fit(train.values.tolist(),label)

    tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
            feature_names = train.columns,
            class_names = ['additive','cooperative'],
            rounded = True, proportion = False,
            precision = 2, filled = True)
    subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])
