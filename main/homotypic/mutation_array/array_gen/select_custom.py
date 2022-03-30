import os
os.chdir("../../../..")
import pandas as pd

import numpy as np
import chip2probe.modeler.mutation as mut
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.pyplot as plt
import chip2probe.training_gen.arranalysis as arr

def get_shifted_pred(df, printlog=False):
    """
    Get prediction that is shifted from coop->add vice versa
    """
    df_sorted = df.sort_values(by=["seqid","main_prob"])
    grouped = df_sorted.groupby('seqid')
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

def make_subselections(df, n):
    """
    We do selection by taking min, max, and take the n middle entries

    The data frame needs to be filtered out with the top candidate before entering
    the function.

    """
    # selection group
    wtdf = df.loc[df["comment"] == "wt"]

    selectg = df.loc[df["comment"] != "wt"] \
        .sort_values(by=["seqid","main_prob"]) \
        .groupby("seqid",as_index=False)
    h_df = pd.DataFrame(selectg.head(1))
    t_df = pd.DataFrame(selectg.tail(1))
    # remove min max rows from each group
    selectg = selectg.apply(lambda g: g.iloc[1:-1]) \
        .groupby("seqid",as_index=False)

    mids = selectg \
        .apply(lambda x: x.sample(n) if n < x['sequence'].count() else x.sample(x['sequence'].count()))

    # we always have wt, head, and tail
    selected = pd.concat([wtdf, h_df, mids, t_df]).sort_values(by=["seqid"]).reset_index(drop=True)
    return selected

def pick_largest_probdif(df, m, n, predcol="main_pred"):
    """
    Get the largest probability difference each group

    Args:
        df: input data frame
        m: how many groups to take
        n: how many elements aside from min max to take from each group
    """
    shiftedpred = get_shifted_pred(df)
    gshifted = shiftedpred.groupby('seqid')["main_prob"] \
            .agg(['max','min']) \
            .assign(diff = lambda d: d["max"]-d["min"]) \
            .sort_values(by=["diff"], ascending=False) \
            .head(m)[['diff']]  # how many groups to take
    selected = df.join(gshifted, on="seqid", how='inner') \
            .drop(["diff"],axis=1) \
            .sort_values(by=["seqid","main_prob"])
    s = make_subselections(selected, n)
    s["select"] = "largest_probdif"
    return s

def pick_consecutive_preddif(df, m, n, predcol="main_pred"):
    grouped = df.groupby('seqid')
    shiftedpred = get_shifted_pred(df, printlog=False)
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
            .head(m) \
            .set_index("seqid")
    selected = df.join(gshift_ct, on="seqid", how='inner') \
                     .drop(["add_count", "coop_count"],axis=1)
    s = make_subselections(selected, n)
    s["select"] = "consecutive_preddif"
    return s

def pick_incos_pred(df, m, n, col1, col2):
    """
    Pick incosistent predictions

    Pick row with different predictions between different prediction columns.
    The priority will be based on the rows that disagree the most. Furthermore,
    we only take rows in the group that disagree.

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
    # here just to get the group with the most inconsistent predictions
    incos = df.loc[df[pred1] != df[pred2]] \
              .assign(diff = lambda x : abs(x[prob1] - x[prob2])) \
              .sort_values(by=["diff"], ascending=False)[["seqid"]] \
              .drop_duplicates() \
              .head(m)
    # get all inconsistent prdictions within group, this might seem redundant
    # since we did the filter above, but I can't think of better way for now
    selected = df.loc[df[pred1] != df[pred2]] \
            .merge(incos, on="seqid", how='inner')
    s = make_subselections(selected, n)
    s["select"] = "incos_pred"
    return s

def pick_similar_affdif(df, n, colname, sim=0.01):
    """
    Pick similiar affinity difference that happen on different side (i.e.
    left and right). All entries are grouped based on this criteria and then
    group selection is made at random.

    Args:
        sim: how close should the number be to be considered similar
    """
    chose = df.copy()
    chose = chose \
        .sort_values(by=["seqid",colname], ascending=True)
    # after
    c = chose.copy()
    c[["comment_t", "%s_t"%colname]] =  chose.groupby('seqid')[["comment",colname]] \
            .shift(-1)
    c["comments"] = c["comment"] + c["comment_t"]
    c = c.loc \
            [
                lambda x: (x["%s_t"%colname] - x["%s"%colname] < sim) &
                 (c["comments"].str.contains("left") & c["comments"].str.contains("right"))
            ][["seqid", "comment", "comment_t"]] \
        .groupby("seqid",as_index=False) \
        .apply(lambda x: x.sample(n) if n < x['seqid'].count() else x.sample(x['seqid'].count())) \
        .reset_index(drop=True)
    cur_c = c[["seqid","comment"]]
    next_c = c[["seqid","comment_t"]].rename(columns = {"comment_t":"comment"})
    toselect = pd.concat([cur_c, next_c]).drop_duplicates()
    chose = chose.merge(toselect, how="inner", on=["seqid","comment"])
    return chose

# affinity
def pick_largest_affdif(df, m, n , sitetarget):
    """
    Args:
        sitetarget: site1 or site2
        m: top m with larger affinity difference to take
        n: how many elements from each group to take, aside from min max
    """

    affdf = pd.DataFrame(df.loc[df["comment"] != "wt"])
    coltarget = "%s_affinity" % sitetarget
    affdf["coltarget"] = sitetarget  #affdf.apply(lambda row: row["%s_affinity" % coltarget],axis=1)
    # only target mutation that affects site1 or site 2
    filt = affdf.apply(lambda row: row["coltarget"] in row["comment"], axis=1)
    affdf = affdf[filt]
    affg = affdf.groupby('seqid')[coltarget] \
            .agg(['max','min']) \
            .assign(diff = lambda d: d["max"]-d["min"]) \
            .sort_values(by=["diff"], ascending=False) \
            .head(m)[["diff"]]
    subselect = affdf.merge(affg, on="seqid", how='inner') \
            .drop(["diff", "coltarget"],axis=1) \
            .sort_values(by=["seqid",coltarget])

    # pick highest, lowest, and some sequences in the middle
    g = subselect.groupby('seqid')
    gh = pd.DataFrame(g.head(1))
    gh["select"] = "largest_affdiff_%s" % sitetarget
    gt = pd.DataFrame(g.tail(1))
    gt["select"] = "largest_affdiff_%s" % sitetarget

    g = g.apply(lambda g: g.iloc[1:-1]) # remove the first (tail -1) and last (head -1) row

    mid = pick_similar_affdif(subselect, n, coltarget)
    mid["select"] = "largest_affdiff_%s" % sitetarget

    # save the wt
    wt_df = df.loc[df["comment"] == "wt"] \
        .merge(subselect[["seqid"]].drop_duplicates(), on="seqid") \
        .assign(select = "largest_affdiff_%s" % sitetarget)
    selected = pd.concat([wt_df, gh,gt, mid]) \
            .sort_values(by=["seqid",coltarget]) \
            .reset_index(drop=True)
    return selected

def pick_positive_ctrl(df, m, n, col1, col2, mincoop=0.5 , coopthres = 0.6, addthres = 0.3):
    """
    Get positive control: high probability and agreement between predictors

    First we look for sequence where wt agrees with high confidence, then
    we take mutants that also agree with high confidence.

    Args:
        minthres: minimum threshold for both predictors, cooperative should be
          above minthres and additive should be below 1-minthres.
        mincoop: how many cooperative at least in the final table, calculted
          as a portion from m

    TODO: see the threshold distribution, make the coop & add more balance
    """
    if coopthres < 0.5 or addthres > 0.5:
        raise ValueError("coopthres needs > 0.5 and addthres < 0.5")

    pred1, pred2 = "%s_pred" % col1, "%s_pred" % col2
    prob1, prob2 = "%s_prob" % col1, "%s_prob" % col2
    # we only use the matching wt
    matchingwt = df.loc[(df["comment"] == "wt") & (df[pred1] == df[pred2])] \
        .loc[
            ((df[pred1] == 1) & (df[prob1] > coopthres) & (df[prob2] > coopthres)) |
            ((df[pred1] == 0) & (df[prob1] < addthres) & (df[prob2] < addthres))
            ]
    coop_m = int(math.ceil(mincoop * m))
    coopsamples = matchingwt.loc[matchingwt["main_pred"] == 1].sample(coop_m)[["seqid"]]
    addsamples = matchingwt.loc[matchingwt["main_pred"] == 0].sample(m - coop_m)[["seqid"]]
    samples = pd.concat([coopsamples,addsamples])
    filtered = df.merge(samples, on=["seqid"], how='inner')
    # do the same filter but among all the filtered groups
    filtered = filtered.loc[
            ((filtered[pred1] == 1) & (filtered[prob1] > coopthres) & (filtered[prob2] > coopthres)) |
            ((filtered[pred1] == 0) & (filtered[prob1] < addthres) & (filtered[prob2] < addthres))
            ]
    filtwt = filtered.loc[filtered["comment"] == "wt"]

    # sample within group
    filtmt = filtered.loc[filtered["comment"] != "wt"] \
        .groupby("seqid") \
        .apply(lambda x: x.sample(n) if n < x['seqid'].count() else x.sample(x['seqid'].count()))

    selected = pd.concat([filtwt, filtmt]) \
        .sort_values("seqid") \
        .reset_index(drop=True)
    selected["select"] = "positive_ctrl"
    return selected

def plot_coop(df, lb,suffix):
    dfsub = df[df["comment"] == "wt"][["Name"]].drop_duplicates()
    ax = plt.axes()
    arr.plot_classified_labels(lb, col1="indiv_median", col2="two_median", log=True, plotnonsignif=False,
                       xlab="M1-M3+M2-M3", ylab="WT-M3", path="labeled_log_one_both.png", title="Cooperative vs independent binding of Ets1-Ets1",
                       labelnames=["cooperative","independent","anticooperative"], axes=ax)
    wtavail = lb.merge(dfsub)
    ax.scatter(np.log(wtavail["indiv_median"]), np.log(wtavail["two_median"]), color="cyan", s=1, label="wt_selected")
    ax.legend()
    plt.savefig("in_%s.png" % suffix)
    plt.clf()

def filter_by_delta(df, lb, ncoop=300,nadd=300):
    """
    get point farthest from the diagonal
    """
    dfnm = df[["Name"]].drop_duplicates()
    wtavail = lb.merge(dfnm)
    wtavail["delta"] = abs(wtavail["two_median"] - wtavail["indiv_median"]) / np.sqrt(2)
    wtcoop = wtavail[wtavail["label"] == "cooperative"].nlargest(ncoop,'delta')[["Name"]]
    wtadd = wtavail[wtavail["label"] == "independent"].nsmallest(nadd,'delta')[["Name"]]
    selected = pd.concat([wtcoop,wtadd])
    seqids_sel = df.merge(selected)[["seqid"]].drop_duplicates()
    return df.merge(seqids_sel, on="seqid")

if __name__ == "__main__":
    df = pd.read_csv("custom_withpred.csv").rename(columns={"name":"Name","Sequence":"sequence"})

    lbled = pd.read_csv("main_nar/output/Ets1Ets1/label_pr/lbled_o1_selected.csv")

    plot_coop(df,lbled,"all")
    # TODO: remove duplicates

    # First we select only groups where predictions and wtlabel are the same
    match_wt = df.loc[(df["comment"] == "wt") & (df["wtlabel"] == df["main_pred"])][["seqid"]] \
        .drop_duplicates()
    mismatch_wt = df.loc[(df["comment"] == "wt") & (df["wtlabel"] != df["main_pred"])][["seqid"]] \
        .drop_duplicates()
    print("Number of correct predicted wt groups %d" % match_wt.shape[0])
    print("Number of incorrect predicted wt groups %d" % mismatch_wt.shape[0])

    pd.set_option("display.max_columns",None)
    # filter the data frame to take only where wt matches
    df = df.merge(match_wt, on="seqid")
    df = filter_by_delta(df,lbled)
    print("wt filtered",df.loc[df["comment"] == "wt"][["sequence"]].drop_duplicates().shape[0])

    muttypes = {"distance": {"ascending":False, "col":"distance"},
                "affinity": {},
                "orientation": {}}
    allpicks = pd.DataFrame()

    for mty in muttypes: #muttypes:
        cur_df = df.loc[df["muttype"] == mty]
        sortby = ["seqid",muttypes[mty]["col"]] if muttypes[mty] else ["seqid"]
        sortasc = [True,muttypes[mty]["ascending"]] if muttypes[mty] else True
        posctrl = pick_positive_ctrl(cur_df, 70, 10, "shape", "main") \
                    .sort_values(by=sortby, ascending=sortasc)
        p1 = pick_largest_probdif(cur_df, 140, 10) \
                    .sort_values(by=["seqid","main_prob"])
        p2 = pick_consecutive_preddif(cur_df, 140, 10) \
                    .sort_values(by=["seqid","main_prob"])
        p3 = pick_incos_pred(cur_df, 140, 10, "shape", "main") \
                    .sort_values(by=sortby, ascending=sortasc)
        p = pd.concat([p1, p2, p3, posctrl])
        if mty == "affinity":
            p_s1 = pick_largest_affdif(cur_df, 110, 10, "site1")
            p_s2 = pick_largest_affdif(cur_df, 110, 10, "site2")
            p = pd.concat([p,p_s1,p_s2])
        allpicks = pd.concat([allpicks, p])
    print(len(allpicks["sequence"].unique()))
    allpicks = allpicks[~allpicks['sequence'].str.contains('GGGGG') & ~allpicks['sequence'].str.contains('CCCCC')]
    allpicks.to_csv("custom_probes_selected.csv", index=False, header=True)

    print("Number of unique sequences %d" % len(allpicks["sequence"].unique()))
    print("Total array spots",len(allpicks["sequence"].unique()) * 24)

    # # Analysis part:
    allpicks = pd.read_csv("custom_probes_selected.csv")
    plot_coop(allpicks,lbled,"after")

    #
    # # get the farthest distance from each group
    # dmax = dist_df \
    #     .sort_values(["seqid", "distance"], ascending=[True,False]) \
    #     .groupby('seqid') \
    #     .head(1)
    # print(dmax)