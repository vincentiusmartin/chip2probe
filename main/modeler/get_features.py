'''
Created on Oct 30, 2019

@author: Vincentius Martin, Farica Zhuang

Pick the best model
'''

import pandas as pd
import os, sys
os.chdir("../..")

from chip2probe.modeler.cooptrain import CoopTrain
# TODO: fix after we finish refactoring probefilter, for now just append the path
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/chip2probe/probe_generator/probefilter")
from sitespredict.imads import iMADS
from sitespredict.imadsmodel import iMADSModel

if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    # select only genomic (i.e. non-custom) sequences
    df = df[~df['name'].str.contains("dist|weak")]
    traindf = CoopTrain(df, corelen=4)

    # using custom imads modoel
    imads_paths = ["input/modeler/imads_model/Ets1_w12_GGAA.model", "input/modeler/imads_model/Ets1_w12_GGAT.model"]
    imads_cores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_paths, imads_cores)]
    imads = iMADS(imads_models, 0.19) # 0.2128

    feature_dict = {
            "distance":{"type":"numerical"},
            "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True},
            "affinity": {"imads":imads}
        }
    traindf.get_feature_all(feature_dict)
