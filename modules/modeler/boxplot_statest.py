'''
Created on Oct 30, 2019

@author: Vincentius Martin, Farica Zhuang

Make some plots for analysis
'''

import pandas as pd
import os
os.chdir("../..")

import chip2probe.modeler.plotmodule as pm
from chip2probe.modeler.cooptrain import CoopTrain

if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    # It is recommended to use the training object so we can get feature specific dataframe
    train = CoopTrain(df, corelen=4)

    # make distace stacked bar
    pm.plot_box_categories(train.df, incols=["distance", "site_str_score", "site_wk_score"])

    pm.plot_grouped_label(train.df, incols=["distance", "site_str_score", "site_wk_score"], by=["label"])
