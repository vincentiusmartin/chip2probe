import os
os.chdir("../..")

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl

if __name__ == "__main__":
    trainingpath = "output/heterotypic/EtsRunx_v1/training_pwm.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    ct = CoopTrain(df)
    pd.set_option("display.max_columns",None)

    rf_param_grid = {
        'n_estimators': [500],
        'max_depth':[5],
        "min_samples_leaf" : [10],
        "min_samples_split" :[10]
    }

    x = ct.get_training_df({"orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}})
    print(x)
    """
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
        "distance,strength":
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
    """
