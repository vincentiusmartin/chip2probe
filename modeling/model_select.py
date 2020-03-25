import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL
from trainingdata.training import Training

from sklearn import ensemble

import pandas as pd

def get_numeric_label(training):
    # hard coded but change add to anti coop / additive when needed
    train = training['label'].map({'cooperative': 1, 'additive': 0})
    return train

def plot_auc(x_train_dict, df, plotname="auc.png"):
    tpr_dict = {key:[] for key in x_train_dict}
    auc_dict = {key:[] for key in x_train_dict}
    acc_dict = {key:[] for key in x_train_dict}

if __name__ == "__main__":
    trainingpath = "train1.tsv"

    df = pd.read_csv(trainingpath, sep="\t")

    t = Training(df, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])

    xtr = t.get_feature_all({
                            "distance":{"type":"numerical", "include":"T"},
                            "orientation":{"positive_cores":["GGAA","GGAT"]},
                            "sitepref":{}
                            })

    x_train = pd.DataFrame(xtr).values.tolist()
    y_train = get_numeric_label(t.df).values
    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=10,random_state=0)
