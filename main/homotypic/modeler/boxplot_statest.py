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
    df = pd.read_csv(trainingpath) #, sep="\t")
    print(df.columns)
    df.rename(columns={'site_str_score': 'Stronger site strength', 'site_wk_score': 'Weaker site strength'}, inplace=True)
    # It is recommended to use the training object so we can get feature specific dataframe
    train = CoopTrain(df, corelen=4)

    # make distace stacked bar
    pl.plot_box_categories(train.df, incols=["distance", "Stronger site strength", "Weaker site strength"], alternative="smaller")

    # pl.plot_grouped_label(train.df, incols=["distance", "site_str_score", "site_wk_score"], by=["label"])
