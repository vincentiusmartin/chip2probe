import os
os.chdir("../..")
import pandas as pd

import chip2probe.modeler.mutation as mut
import matplotlib.pyplot as plt

def get_shifted_pred(df, printlog=True):
    grouped = df.groupby('seqid')
    shiftcount = 0
    wrongpred_count = 0
    shiftpred = []
    for name, group in grouped:
        first = group.head(1).iloc[0]
        # first, ignore if prediction is different with wt
        if first["wtlabel"] != first["main_pred"]:
            wrongpred_count += 1
            continue
        last = group.tail(1).iloc[0]
        if first["main_pred"] != last["main_pred"]:
            shiftpred.extend(group.to_dict('records'))
            shiftcount += 1
    if printlog:
        print("Number of original wt sequences %d" % len(grouped))
        print("Number of wrong wt predictions %d" % wrongpred_count)
        print("Number of changed predictions:  %d" % shiftcount)
    return pd.DataFrame(shiftpred)

def pick_largest_probdif(df, n, predcol="main_pred"):
    shiftedpred = get_shifted_pred(df)
    gshifted = shiftedpred.groupby('seqid')["main_prob"] \
            .agg(['max','min']) \
            .assign(diff = lambda d: d["max"]-d["min"]) \
            .sort_values(by=["diff"], ascending=False) \
            .head(n)[["diff"]]  # how many groups to take
    selected = df.join(gshifted, on="seqid", how='inner').drop(["diff"],axis=1)
    selected["select"] = "largest_probdif"
    return selected

def pick_consecutive_preddif(df, n, predcol="main_pred"):
    grouped = df.groupby('seqid')
    shiftedpred = get_shifted_pred(df, printlog=True)
    gshifted = shiftedpred.groupby(["seqid",predcol])["sequence"] \
            .agg(["count"]).reset_index()
    gshift_addct = gshifted.loc[gshifted[predcol] == 0] \
            .rename(columns={"count":"add_count"}).drop([predcol],axis=1) \
            .set_index("seqid")
    gshift_coopct = gshifted.loc[gshifted[predcol] == 1] \
            .rename(columns={"count":"coop_count"}).drop([predcol],axis=1) \
            .set_index("seqid")
    gshift_ct = gshift_addct.join(gshift_coopct, on="seqid", how='inner') \
            .reset_index() \
            .sort_values(by=["coop_count"], ascending=False) \
            .head(n) \
            .set_index("seqid")
    selected = df.join(gshift_ct, on="seqid", how='inner') \
                     .drop(["add_count", "coop_count"],axis=1)
    selected["select"] = "consecutive_preddif"
    return selected

def pick_incos_pred(df, n, col1, col2):
    """
    Pick incosistent predictions

    Pick row with different predictions between different prediction columns.
    The priority will be based on the rows that disagree the most.

    Args:
        df: data frame input
        n: how many top entries to take
        predcols: list of 2 containing different predictions where the prediction
            will have "pred" as suffix and "prob" suffix for probability.
    Returns:
        selected: a data frame of selected rows
    """
    pred1, pred2 = "%s_pred" % col1, "%s_pred" % col2
    prob1, prob2 = "%s_prob" % col1, "%s_prob" % col2
    incos = df.loc[df[pred1] != df[pred2]] \
              .assign(diff = lambda x : abs(x[prob1] - x[prob2])) \
              .sort_values(by=["diff"], ascending=False) \
              .drop_duplicates() \
              .head(n)[["seqid"]] \
              .set_index("seqid")
    selected = df.join(incos, on="seqid", how='inner')
    selected["select"] = "incos_pred"
    selected.to_csv("blaaa.csv", index=False, header=True)
    return selected

if __name__ == "__main__":
    df = pd.read_csv("output/custom_sequences/custom_distance_withpred.csv")
    print(df.shape)
    print((df.drop_duplicates()).shape)
    plt.scatter(df['main_prob'], df['shape_prob'],s=1, c='blue')
    plt.savefig("corr.eps")
    r = df['main_prob'].corr(df['shape_prob'])
    print(r*r)
    # TODO: remove duplicates


    #p1 = pick_largest_probdif(df, 20)
    # p2 = pick_consecutive_preddif(df, 20)
    # p3 = pick_incos_pred(df, 20, "shape","main")
    # p = pd.concat([p1,p2,p3])
    # print(p.shape)
    # p.to_csv("blaaa.csv", index=False, header=True)
