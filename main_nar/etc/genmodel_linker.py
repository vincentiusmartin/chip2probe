import os
os.chdir("..")

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl
import pickle
from sklearn import ensemble, tree
import subprocess

import chip2probe.modeler.features.sequence as seq


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

    s = seq.Sequence(ct.df, {"seqin":2, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":"linker"})
    s.get_feature()

    rf_param_grid = {
        'n_estimators': [1000], #[500,750,1000],
        'max_depth': [10], #[5,10,15],
        "min_samples_leaf": [5], #[5,10,15],
        "min_samples_split" : [5] #[5,10,15]
    }

    best_models = {
        "2":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "sequence":{"seqin":2, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":"linker"}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "3":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "sequence":{"seqin":3, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":"linker"}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "4":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "sequence":{"seqin":4, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":"linker"}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
        "6":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "sequence":{"seqin":4, "poscols":['%s_pos'%s1,'%s_pos'%s2], "namecol":"Name", "smode":"linker"}
                }, label_map={'cooperative': 1, 'independent': 0})
            ).run_all(),
    }

    pl.plot_model_metrics(best_models, path="auc_separate_str.png" % basepath, cvfold=10, score_type="auc", varyline=True, title="Average ROC Curves for Ets1-Runx1")

    # tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
    #         feature_names = train.columns,
    #         class_names = ['additive','cooperative'],
    #         rounded = True, proportion = False,
    #         precision = 2, filled = True)
    # subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])
