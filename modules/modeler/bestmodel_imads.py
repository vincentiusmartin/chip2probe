'''
Created on Oct 30, 2019

@author: Vincentius Martin, Farica Zhuang

Pick the best model
'''

import pandas as pd
import os, sys
os.chdir("../..")

import chip2probe.modeler.plotlib as pl
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
# TODO: fix after we finish refactoring probefilter, for now just append the path
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/chip2probe/probe_generator/probefilter")
from sitespredict.imads import iMADS
from sitespredict.imadsmodel import iMADSModel

if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    # select only genomic (i.e. non-custom) sequences
    df = df[~df['name'].str.contains("dist|weak")]
    cooptr = CoopTrain(df, corelen=4)

    rf_param_grid = {
        'n_estimators': [500, 1000, 1500],
        'max_depth':[5, 10, 15],
        "min_samples_leaf" : [10, 15, 20],
        "min_samples_split" :[10, 15 ,20]
    }

    # using custom imads model
    imads8_paths = ["input/modeler/imads_model/Ets1_w8_GGAA.model", "input/modeler/imads_model/Ets1_w8_GGAT.model"]
    imads8_cores = ["GGAA", "GGAT"]
    imads8_models = [iMADSModel(path, core, 8, [1, 2, 3]) for path, core in zip(imads8_paths, imads8_cores)]
    imads8 = iMADS(imads8_models, 0.19) # 0.2128

    imads12_paths = ["input/modeler/imads_model/Ets1_w12_GGAA.model", "input/modeler/imads_model/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19) # 0.2128

    best_models = {
        "imads8":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=cooptr.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"imads": imads8},
                    "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                })
            ).run_all(),
        "imads12":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=cooptr.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"imads": imads12},
                    "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                })
            ).run_all(),
        "imads20":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=cooptr.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {},
                    "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                })
            ).run_all()
    }
    # import pickle
    # best_models = pickle.load(open( "bm.pickle", "rb" ) )
    pl.plot_model_metrics(best_models, cvfold=10, score_type="auc")
