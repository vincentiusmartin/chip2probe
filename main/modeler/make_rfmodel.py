'''
Created on Oct 30, 2019

@author: Vincentius Martin, Farica Zhuang

Pick the best model
'''

import pandas as pd
import os, sys
import pickle
from sklearn import ensemble, tree
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
    ct = CoopTrain(df, corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"])

    # using custom imads model
    imads_paths = ["input/modeler/imads_model/Ets1_w12_GGAA.model", "input/modeler/imads_model/Ets1_w12_GGAT.model"]
    imads_cores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_paths, imads_cores)]
    imads = iMADS(imads_models, 0.19) # 0.2128

    # get the features from the CoopTrain class
    feature_dict = {
            "distance":{"type":"numerical"},
            "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True},
            "affinity": {"imads":imads}
        }
    train = ct.get_feature_all(feature_dict, aslist=True)
    label = ct.get_numeric_label({'cooperative': 1, 'additive': 0})

    tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
            feature_names = xt_df.columns,
            class_names = ['additive','cooperative'],
            rounded = True, proportion = False,
            precision = 2, filled = True)
    subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])
