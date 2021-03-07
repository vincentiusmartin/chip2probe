import os
os.chdir("../../../..")

import chip2probe.training_gen.arranalysis as arr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_subrows(path, inlist):
    df = pd.read_csv(path)
    df = df[df["Name"].isin(inlist)]
    return df


if __name__ == "__main__":
    pd.set_option("display.max_columns", None)

    indiv = get_subrows("output/homotypic/custom_sequences/20nM/indiv.csv", [2642,2647])
    two = get_subrows("output/homotypic/custom_sequences/20nM/two.csv", [2642,2647])
    arr.plot_ori_inconsistency(indiv, two, namecol="Name", prefix_path="all", log=True, allplot=False)
