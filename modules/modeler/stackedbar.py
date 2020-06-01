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
    train = CoopTrain(df, corelen=4)

    # make distace stacked bar
    pm.stacked_bar_categories(train.df, "distance", path="dist_stackedbar.png")

    # get stacked bar of ratio between different distance
    pm.stacked_bar_categories(train.df, "distance", path="dist_stackedbar_avg.png", avg=True)
