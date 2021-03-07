'''
Created on Oct 30, 2019

@author: Vincentius Martin, Farica Zhuang

Make some plots for analysis
'''

import pandas as pd
import os
os.chdir("../../..")

import chip2probe.modeler.plotlib as pl
from chip2probe.modeler.cooptrain import CoopTrain

if __name__ == "__main__":
    trainingpath = "output/homotypic/training/training_pwm.csv"
    df = pd.read_csv(trainingpath) # , sep="\t")
    df['label'] = df['label'].replace({"additive":"independent"})
    train = CoopTrain(df, corelen=4)

    # make distace stacked bar
    pl.plot_stacked_categories(train.df, "orientation", path="dist_stackedbar.png")

    # get stacked bar of ratio between different distance
    pl.plot_stacked_categories(train.df, "distance", path="distance_ets_ets.png", ratio=True)
