import pandas as pd

from sitespredload import imads_ets, imads_runx, kp_ets, kp_runx
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt

pd.set_option("display.max_columns",None)
if __name__ == "__main__":
    df = pd.read_csv("output/mutall_etsets_with_imads.csv")
    s1deltcut, s2deltcut = 0.1, 0.1
    dfgrp = df.groupby("id")

    allgrps = []
    for i,g in dfgrp:
        wt_g = g[g["comment"] == "wt"]
        if wt_g.empty:
            print("empty",i)
            continue
        wt_g = wt_g.iloc[0]
        for s in ["s1","s2"]:
            g["%s_score_delta"%s] = g["%s_score"%s] - wt_g["%s_score"%s]
            g["%s_imads_delta"%s] = g["%s_imads"%s] - wt_g["%s_imads"%s]
        allgrps.append(g)
    resdf = pd.concat(allgrps)

    alldelta = {}
    for s in ["s1","s2"]:
        x, y = resdf["%s_score_delta"%s].tolist(), resdf["%s_imads_delta"%s].tolist()
        x = [[elm] for elm in x]
        model = LinearRegression()
        model.fit(x, y)
        alldelta[s] = abs(y - model.predict(x))

    resdf["s1_ldist"] = alldelta["s1"]
    resdf["s2_ldist"] = alldelta["s2"]

    resdf = resdf[
            (resdf["comment"] == "wt") |
                (
                    (np.sign(resdf["s1_score"]) == np.sign(resdf["s1_imads"])) &
                    (np.sign(resdf["s2_score"]) == np.sign(resdf["s2_imads"])) &
                    (resdf["s1_ldist"] < s1deltcut) & (resdf["s2_ldist"] < s2deltcut)
                )
            ]
    resdf.to_csv("mutall_etsets_delta.csv")
    print(resdf.shape[0])

    resdf[["s1_score","s2_score"]].boxplot()
    plt.savefig("pwmbox.png")
    plt.clf()

    for s in ["s1","s2"]:
        print("Rsq", resdf["%s_score_delta"%s].corr(resdf["%s_imads_delta"%s])**2)
        x, y = resdf["%s_score_delta" % s].tolist(), resdf["%s_imads_delta" % s].tolist()
        plt.scatter(x,y,color="blue",s=0.8)
        slope, intercept = np.polyfit(x, y, 1)
        abline_values = [slope * i + intercept for i in x]
        plt.xlabel("PWM delta")
        plt.ylabel("iMADS delta")
        plt.plot(x, abline_values, color='red', linestyle='dashed')
        plt.title("Delta %s" % s)
        plt.savefig('scatter_%s.png'%s)
        plt.clf()
