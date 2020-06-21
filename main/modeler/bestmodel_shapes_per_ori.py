'''
Created on June 9, 2020

@author: Vincentius Martin, Farica Zhuang

Pick the best shape model for different orientation
'''

import pandas as pd
import os
os.chdir("../..")

import chip2probe.modeler.plotlib as pl
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
from chip2probe.modeler.dnashape import DNAShape

if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    # select only genomic (i.e. non-custom) sequences
    df = df[~df['name'].str.contains("dist|weak")]
    cooptr = CoopTrain(df, corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"])
    #
    x_ori = cooptr.get_feature("orientation", {"positive_cores":["GGAA", "GGAT"]})
    df["orientation"] = pd.DataFrame(x_ori)["ori"]

    rf_param_grid = {
        'n_estimators': [500],
        'max_depth':[5],
        "min_samples_leaf" : [10],
        "min_samples_split" :[10]
    }

    orientations = ["HH", "TT", "HT/TH"]
    ori = "HH"


    df_ht = dftrain #[dftrain["orientation"] == "HT/TH"] #[dftrain["distance"] % 2 == 0] #[dftrain["orientation"] == "HT/TH"]
    dftr.to_csv("train1.tsv",sep="\t")

    for ori in orientations:
        print(ori)

    # # TODO: choose per orientation
    # best_models = {"top10":
    #           BestModel(clf="sklearn.ensemble.RandomForestClassifier",
    #                       param_grid=rf_param_grid,
    #                       train_data=cooptr.get_training_df({
    #                               "distance":{"type":"numerical"},
    #                               "shape": {"seqin":2, "smode":"positional", "direction":"inout"},
    #                               "shape": {"ds":ds, "seqin":-3, "smode":"strength"},
    #                               #"orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    #                           })
    #             )#.run_all()
    # }
    # pl.plot_model_metrics(best_models, cvfold=10, score_type="auc")
