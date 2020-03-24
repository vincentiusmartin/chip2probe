import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL
from trainingdata.training import Training

import pandas as pd



if __name__ == "__main__":
    trainingpath = "train1.tsv"

    df = pd.read_csv(trainingpath, sep="\t")

    t = Training(df, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])
