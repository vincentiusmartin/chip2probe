'''
Created on Oct 30, 2019

@author: Vincentius Martin, Farica Zhuang

Pick the best model
'''

import pandas as pd
import os, sys
import pickle
from sklearn import ensemble, tree
import chip2probe.modeler.plotlib as pl
os.chdir("../../..")

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
# TODO: fix after we finish refactoring probefilter, for now just append the path
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

if __name__ == "__main__":
    trainingpath =  "input/modeler/training_data/training_p01_adjusted_ets1.tsv"
    df = pd.read_csv(trainingpath,sep="\t")
    # select only genomic (i.e. non-custom) sequences
    # df = df[~df['name'].str.contains("dist|weak")]
    ct = CoopTrain(df, corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"])

    # using custom imads model
    imads_paths = ["input/site_models/imads_models/Ets1_w12_GGAA.model", "input/site_models/imads_models/Ets1_w12_GGAT.model"]
    imads_cores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_paths, imads_cores)]
    imads = iMADS(imads_models, 0.19) # 0.2128

    # get the features from the CoopTrain class
    feature_dict = {
            "distance":{"type":"numerical"},
            "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True},
            "affinity": {"imads":imads}
        }
    train = ct.get_feature_all(feature_dict)
    label = ct.get_numeric_label({'cooperative': 1, 'additive': 0})

    rf_param_grid = {
        'n_estimators': [1000], #[500,750,1000],
        'max_depth': [10], #[5,10,15],
        "min_samples_leaf": [10], #[5,10,15],
        "min_samples_split" : [10]#[5,10,15]
    }


    best_models = {
        "all":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df(feature_dict)
            ).run_all()
    }
    pl.plot_model_metrics(best_models, cvfold=10, score_type="auc", varyline=True, title="Average ROC Curves for Runx1-Ets1")

    # # This is usually made based on the best model
    rf = ensemble.RandomForestClassifier(n_estimators=1000, max_depth=10, min_samples_leaf=10, min_samples_split=10)
    rf.fit(train,label)
    model_name = "dist_ori_12merimads.sav"
    pickle.dump(rf, open(model_name, 'wb'))
    print("Model saved in %s" % model_name)

    # tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
    #         feature_names = xt_df.columns,
    #         class_names = ['additive','cooperative'],
    #         rounded = True, proportion = False,
    #         precision = 2, filled = True)
    # subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])
