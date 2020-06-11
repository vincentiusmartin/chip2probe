'''
Created on June 9, 2020

@author: Vincentius Martin, Farica Zhuang

Pick the best model
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
    cooptr = CoopTrain(df, corelen=4)

    shapepath = "input/modeler/dnashape/0"
    ds = DNAShape(shapepath)

    rf_param_grid = {
        'n_estimators': [500, 1000, 1500],
        'max_depth':[5, 10, 15],
        "min_samples_leaf" : [10, 15, 20],
        "min_samples_split" :[10, 15 ,20]
    }

    best_models = {"top10":
              BestModel(clf="RF",
                          param_dict=rf_param_grid,
                          train_data=cooptr.get_training_df({
                                  #"distance":{"type":"numerical"},
                                  "shape": {"dnashape":ds, "seqin":5, "smode":"strength"}
                                  #"shape": {"ds":ds, "seqin":-3, "smode":"strength"}
                                  #"orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                              }),
                           topn=10
                ).run_all()
    }
