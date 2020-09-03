import pandas as pd
import itertools
import os

from inputdict import param
import imads_train as mt

import libsvm.svmutil as svmutil


if __name__ == "__main__":
    num_workers = os.cpu_count()

    pd.set_option('display.max_columns', None)
    cores_centered = mt.init_train_matrix(param)

    mt.genmodel_gridsearch(cores_centered, param["grid"], param["numfold"], param["kmers"], num_workers,
                           logit=param["logit"], tfname=param["tfname"], modelwidth=param["width"],
                           outdir=param["outdir"])
