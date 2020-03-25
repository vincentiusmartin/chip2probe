import sys
#sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL
sys.path.append("/Users/faricazjj/Desktop/homotf/chip2probe")
from trainingdata.training import Training
from best_model import BestModel

from sklearn import ensemble

import pandas as pd

def plot_auc(x_train_dict, df, plotname="auc.png"):
    tpr_dict = {key:[] for key in x_train_dict}
    auc_dict = {key:[] for key in x_train_dict}
    acc_dict = {key:[] for key in x_train_dict}

if __name__ == "__main__":
    #trainingpath = "train1.tsv"
    trainingpath = "trainingdata/training_new.csv"

    df = pd.read_csv(trainingpath, sep=",")

    t = Training(df, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])

    train_data = t.get_training_df({
                            "distance":{"type":"numerical", "include":"T"},
                            "orientation":{"positive_cores":["GGAA","GGAT"]},
                            "sitepref":{}
                            })
    param_dict = {
    				'n_estimators': [i for i in range(2,21)],
    				'max_depth': [i for i in range(100,2001,100)]
    			}
    rf = BestModel(clf="RF", param_dict=param_dict, topn=10, train_data=train_data).run_all()
