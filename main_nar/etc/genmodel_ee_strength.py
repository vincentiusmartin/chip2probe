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
    basepath = "output/Ets1Ets1"
    trainingpath = "output/Ets1Ets1/training/train_ets1_ets1.tsv"

    df = pd.read_csv(trainingpath, sep="\t")
    df["avg_str"] = (df["site_wk_score"] + df["site_str_score"]) / 2
    ct = CoopTrain(df)
    pd.set_option("display.max_columns",None)


    rf_param_grid = {
        'n_estimators': [750,1000], #[500,750,1000],
        'max_depth': [5,10], #[5,10,15],
        "min_samples_leaf": [5,10,15], #[5,10,15],
        "min_samples_split" : [5,10,15] #[5,10,15]
    }

    best_models = {
        "Weaker site strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "affinity": {"colnames": ["site_str_score"]}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "Stronger site strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "affinity": {"colnames": ["site_wk_score"]}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,wkr_strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ["site_str_score"]},
                    "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation,str_strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ["site_wk_score"]},
                    "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
            "distance,orientation,avg_strength":
                BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                  param_grid=rf_param_grid,
                  train_data=ct.get_training_df({
                        "distance":{"type":"numerical"},
                        "affinity": {"colnames": ["avg_str"]},
                        "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
                    }, label_map={'cooperative': 1, 'independent': 0})
                ).run_all(),
            "distance,orientation,strength":
                BestModel(clf="sklearn.ensemble.RandomForestClassifier",
                  param_grid=rf_param_grid,
                  train_data=ct.get_training_df({
                        "distance":{"type":"numerical"},
                        "affinity": {"colnames": ("site_str_score","site_wk_score")},
                        "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
                    }, label_map={'cooperative': 1, 'independent': 0})
                ).run_all()
    }

    pl.plot_model_metrics(best_models, path="%s/auc_separate_str.pdf" % basepath, cvfold=10, score_type="auc", varyline=True, title="ROC Curves for Ets1-Ets1",interp=True)

    # tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
    #         feature_names = train.columns,
    #         class_names = ['additive','cooperative'],
    #         rounded = True, proportion = False,
    #         precision = 2, filled = True)
    # subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])