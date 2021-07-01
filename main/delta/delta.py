import os
os.chdir("../..")

import pandas as pd
import numpy as np

from chip2probe.modeler.cooptrain import CoopTrain
from sklearn.model_selection import cross_validate
from sklearn.model_selection import train_test_split

from sklearn.ensemble import RandomForestRegressor
import sklearn.metrics as metrics

from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

"""
1. Using distance, affinity, orientation
RandomForestRegressor: small min_samples_leaf seems to be the strongest predictor, for (n_estimators=200, max_depth=30, random_state=0, min_samples_leaf=2, min_samples_split=2) r2_train=0.77
2. Using dist, ori, shape, seq
R2 = 0.78
"""

def svc_param_selection(X, y, nfolds):
    Cs = [0.001, 0.01, 0.1, 1, 10]
    gammas = [0.001, 0.01, 0.1, 1]
    epsilon = [0.0001, 0.001, 0.01, 0.1, 0.5, 1]
    kernel = ['linear', 'rbf']
    param_grid = {'C': Cs, 'gamma' : gammas, 'kernel':kernel, 'epsilon':epsilon}
    grid_search = GridSearchCV(SVR(), param_grid, cv=nfolds)
    grid_search.fit(X, y)
    grid_search.best_params_
    return grid_search.best_params_

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    dfdelta = pd.read_csv("main_nar/output/Ets1Ets1/label_pr/lbled_o1_selected.csv")
    dfdelta["delta"] = dfdelta["two_median"] - dfdelta["indiv_median"]
    dfdelta = dfdelta[["Name","delta"]]

    dft = pd.read_csv("main_nar/output/Ets1Ets1/training/train_ets1_ets1.tsv", sep="\t")
    dft = dft.merge(dfdelta, on="Name")

    ct = CoopTrain(dft)
    feature_dict = {
        "distance":{"type":"numerical"},
        "affinity": {"colnames": ("site_str_score","site_wk_score")}, # strength
        "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}},
        "shape_in":{"seqin":5, "poscols":['site_str_pos','site_wk_pos'], "smode":"relative"},
        "shape_out":{"seqin":-2, "poscols":['site_str_pos','site_wk_pos'], "smode":"relative"},
        "sequence_in":{"seqin":5, "poscols":['site_str_pos','site_wk_pos'], "namecol":"Name", "smode":"relative"},
        "sequence_out":{"seqin":-3, "poscols":['site_str_pos','site_wk_pos'], "namecol":"Name", "smode":"relative"}
    }
    X = ct.get_feature_all(feature_dict).values.tolist()
    print(X)
    # import sys
    # sys.exit()
    X = StandardScaler().fit_transform(X)
    ytrue = ct.df["delta"].values.tolist()

    regr = RandomForestRegressor(n_estimators=200, max_depth=30, random_state=0, min_samples_leaf=2, min_samples_split=2)
    # param = svc_param_selection(X, ytrue, 5)
    # regr = SVR(**param)
    # print("Best",param)

    regr.fit(X,ytrue)
    ypred = regr.predict(X)
    print("Training accuracy",metrics.r2_score(ytrue, ypred))

    cv_results = cross_validate(regr, X, ytrue, cv=10, scoring=('r2', 'neg_mean_squared_error'))
    print(cv_results)
    print(np.mean(cv_results["test_r2"]))
