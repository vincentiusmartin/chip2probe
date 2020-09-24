import os
os.chdir("../..")

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl

if __name__ == "__main__":
    trainingpath = "output/heterotypic/EtsRunx_v1/training_p06.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    ct = CoopTrain(df)

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
        "strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "affinity": {"colnames": ("ets_score","runx_score")}
                })
            ).run_all(),
        "distace,strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("ets_score","runx_score")}
                })
            ).run_all()
    }
    # import pickle
    # best_models = pickle.load(open( "bm.pickle", "rb" ) )
    pl.plot_model_metrics(best_models, cvfold=10, score_type="auc", varyline=True)
