import os
os.chdir("../..")
import pandas as pd

import chip2probe.modeler.mutation as mut
import matplotlib.pyplot as plt

if __name__ == "__main__":
    df = pd.read_csv("output/custom_sequences/custom_distance_withpred.csv")
    plt.scatter(df['main_prob'], df['shape_prob'],s=1, c='blue')
    plt.savefig("corr.eps")
    r = df['main_prob'].corr(df['shape_prob'])
    print(r*r)

    mut.pick_custom_dist(df)
