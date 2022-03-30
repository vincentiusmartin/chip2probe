import os

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl
from sklearn import ensemble, tree
import subprocess
import pickle

pd.set_option("display.max_columns",None)
if __name__ == "__main__":
    basepath = "output/Ets1Ets1_v2"
    trainingpath = "%s/training/train_ets1_ets1.tsv" % basepath
    df = pd.read_csv(trainingpath,sep="\t")
    ct = CoopTrain(df)
    pd.set_option("display.max_columns",None)

    rf_param_grid = {
        'n_estimators': [500] ,#[500,750,1000],
        'max_depth': [10],#[5,10,15],
        "min_samples_leaf": [15],#[5,10,15],
        "min_samples_split" : [10],#[5,10,15]
    }

    # {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}

    best_models = {
        "strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "affinity": {"colnames": ("site_str_score","site_wk_score")}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "orientation":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "distance,orientation":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "strength,distance":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("site_str_score","site_wk_score")}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "strength,distance,orientation":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("site_str_score","site_wk_score")},
                    "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all()
    }

    pl.plot_model_metrics(best_models, path="%s/model/auc_all.png" % basepath, cvfold=10, score_type="auc", varyline=True, title="Average ROC Curves for Ets1-Ets1", interp=True)

    rf = best_models["distance,orientation,strength"][1]

    train = ct.get_feature_all({
        "distance":{"type":"numerical"},
        "affinity": {"colnames": ("site_str_score","site_wk_score")},
        "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
    })
    label = ct.get_numeric_label({'cooperative': 1, 'independent': 0})

    rf.fit(train,label)
    model_name = "%s/model/ets1_ets1_rfmodel.sav" % basepath
    pickle.dump(rf, open(model_name, 'wb'))
    print("Model saved in %s" % model_name)

    """
    feature_dict = {
        "distance":{"type":"numerical"},
        "affinity": {"colnames": ("ets_score","runx_score")},
        "orientation": {"relative":False, "one_hot":True, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}
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
    """
