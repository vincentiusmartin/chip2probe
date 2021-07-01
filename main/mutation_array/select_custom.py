import pandas as pd
import numpy as np
import math
import os

def get_count(df):
    wt_unique = df[df["comment"] == "wt"].drop_duplicates("id")
    coop_wt = wt_unique.loc[wt_unique["wtlabel"] == 1]
    indep_wt =  wt_unique.loc[wt_unique["wtlabel"] == 0]

    for predtype in ["main_pred", "shape_pred"]:
        coop_match = coop_wt.loc[coop_wt["wtlabel"] == coop_wt[predtype]][["id"]] \
            .drop_duplicates()
        indep_match = indep_wt.loc[indep_wt["wtlabel"] == indep_wt[predtype]][["id"]] \
            .drop_duplicates()
        print("Number of correct predicted cooperative wt, %s: %d/%d" % (predtype,coop_match.shape[0],coop_wt.shape[0]))
        print("Number of correct predicted independent wt, %s: %d/%d" % (predtype,indep_match.shape[0],indep_wt.shape[0]))

def pick_positive_ctrl(df, m, n, predcol1, predcol2, mincoop=0.5, coopthres = 0.7, indepthres = 0.3):
    """
    Get positive control: high probability and agreement between predictors

    First we look for sequence where wt agrees with high confidence, then
    we take mutants that also agree with high confidence.

    Args:
        minthres: minimum threshold for both predictors, cooperative should be
          above minthres and additive should be below 1-minthres.
        mincoop: how many cooperative at least in the final table, calculted
          as a portion from m

    """
    if coopthres < 0.5 or indepthres > 0.5:
        raise ValueError("coopthres needs > 0.5 and indepthres < 0.5")

    pred1, pred2 = "%s_pred" % predcol1, "%s_pred" % predcol2
    prob1, prob2 = "%s_proba" % predcol1, "%s_proba" % predcol2

    # we only use the matching wts that fulfill threshold
    matchingwt = df.loc[(df["comment"] == "wt") & (df[pred1] == df[pred2])] \
        .loc[
            ((df[pred1] == 1) & (df[prob1] > coopthres) & (df[prob2] > coopthres)) |
            ((df[pred1] == 0) & (df[prob1] < indepthres) & (df[prob2] < indepthres))
            ]
    coop_m = int(math.ceil(mincoop * m))
    coop_sample = min(matchingwt.loc[matchingwt["main_pred"] == 1].shape[0], coop_m)
    indep_sample = min(matchingwt.loc[matchingwt["main_pred"] == 0].shape[0], m - coop_m)

    coopsamples = matchingwt.loc[matchingwt["main_pred"] == 1].sample(coop_sample)[["id"]]
    indepsamples = matchingwt.loc[matchingwt["main_pred"] == 0].sample(indep_sample)[["id"]]
    samples  = pd.concat([coopsamples, indepsamples])
    filtered = df.merge(samples, on=["id"], how='inner')
    # do the same filter but among all the filtered groups
    filtered = filtered.loc[
            ((filtered[pred1] == 1) & (filtered[prob1] > coopthres) & (filtered[prob2] > coopthres)) |
            ((filtered[pred1] == 0) & (filtered[prob1] < indepthres) & (filtered[prob2] < indepthres))
            ]
    filtwt = filtered.loc[filtered["comment"] == "wt"]

    # sample within group
    filtmt = filtered.loc[filtered["comment"] != "wt"] \
        .groupby("id") \
        .apply(lambda x: x.sample(n) if n < x['id'].count() else x.sample(x['id'].count()))

    selected = pd.concat([filtwt, filtmt]) \
        .sort_values("id") \
        .reset_index(drop=True)
    selected["select"] = "positive_ctrl"

    return selected

def pick_incos_pred(df, m, n, col1="main", col2="shape"):
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
    prob1, prob2 = "%s_proba" % col1, "%s_proba" % col2
    # here just to get the group with the most inconsistent predictions
    incos = df.loc[df[pred1] != df[pred2]] \
              .assign(diff = lambda x : abs(x[prob1] - x[prob2])) \
              .sort_values(by=["diff"], ascending=False)[["id"]] \
              .drop_duplicates() \
              .head(m)
    # get all inconsistent prdictions within group, this might seem redundant
    # since we did the filter above, but I can't think of better way for now
    selected = df.loc[df[pred1] != df[pred2]] \
            .merge(incos, on="id", how='inner')
    s = make_subselections(selected, n)
    s["select"] = "incos_pred"
    return s

def get_shifted_pred(df, mty):
    """
    Get prediction that is shifted from coop->add vice versa
    """
    grouped = df.groupby(by=["id"])
    shiftpred = []
    shiftcount = 0
    for name, group in grouped:
        if group.shape[0] == 1:
            # not enough member
            continue
        flag = False
        if mty == "distance":
            group = group.sort_values("distance")
            if group.head(1).iloc[0]["main_pred"] != group.tail(1).iloc[0]["main_pred"]:
                flag=True
        else:
            wtpred = group[group['comment'] == "wt"].iloc[0]["main_pred"]
            mtpred = group[(group['comment'] != "wt") & (group['main_pred'] != wtpred)]
            if not mtpred.empty:
                # for strength, set the minimum number of mutants
                if (mty == "orientation") or (mty == "strength" and mtpred.shape[0] > 1):
                    flag = True
        if flag:
            shiftpred.extend(group.to_dict('records'))
    return pd.DataFrame(shiftpred)

def make_subselections(df, n, col="main_proba"):
    """
    We do selection by taking min, max, and take the n middle entries

    The data frame needs to be filtered out with the top candidate before entering
    the function.

    """
    # selection group
    wtdf = df.loc[df["comment"] == "wt"]

    selectg = df.loc[df["comment"] != "wt"] \
        .sort_values(by=["id",col]) \
        .groupby("id",as_index=False)
    h_df = pd.DataFrame(selectg.head(1))
    t_df = pd.DataFrame(selectg.tail(1))
    # remove min max rows from each group
    selectg = selectg.apply(lambda g: g.iloc[1:-1]) \
        .groupby("id",as_index=False)

    mids = selectg \
        .apply(lambda x: x.sample(n) if n < x['id'].count() else x.sample(x['id'].count()))

    # we always have wt, head, and tail
    selected = pd.concat([wtdf, h_df, mids, t_df]).sort_values(by=["id"]).reset_index(drop=True)
    return selected

def pick_largest_probdif(df, m, n, muttype, predcol="main_pred"):
    """
    Get the largest probability difference each group

    Args:
        df: input data frame
        m: how many groups to take
        n: how many elements aside from min max to take from each group
    """
    shiftedpred = get_shifted_pred(df, muttype)
    gshifted = shiftedpred.groupby('id')["main_proba"] \
            .agg(['max','min']) \
            .assign(diff = lambda d: d["max"]-d["min"]) \
            .sort_values(by=["diff"], ascending=False) \
            .head(m) \
            .reset_index()[['id']] \
            .drop_duplicates()
    # we need to get more sequences
    if gshifted.shape[0] < m:
        gshifted_ids = set(gshifted["id"])
        notselected = set(df["id"]) - gshifted_ids
        notsel_df = pd.DataFrame({'id':list(notselected)})
        more_rows = df.merge(notsel_df) \
            .groupby('id')["main_proba"] \
            .agg(['max','min']) \
            .assign(diff = lambda d: d["max"]-d["min"]) \
            .sort_values(by=["diff"], ascending=False) \
            .head(m-gshifted.shape[0]) \
            .reset_index()[['id']] \
            .drop_duplicates()
        allids = gshifted_ids.union(set(more_rows['id']))
        gshifted = pd.DataFrame({'id':list(allids)})
    selected = df.merge(gshifted, on="id", how='inner') \
            .sort_values(by=["id","main_proba"])
    s = make_subselections(selected, n)
    s["select"] = "largest_probdif"
    return s

# affinity
def pick_small_affdif(df, m, n ,sitetarget, affdelta = 2):
    """
    Args:
        sitetarget: site1 or site2
        m: top m with larger affinity difference to take
        n: how many elements from each group to take, aside from min max
    """
    coltarget = "%s_score" % sitetarget

    # only target mutation that affects site1 or site 2
    affdf = df[df['comment'].str.contains(sitetarget) | (df["comment"] == "wt")]
    affgrp = affdf.groupby("id") \
                 .filter(lambda g: 0 in g["main_pred"].values and 1 in g["main_pred"].values) \
                 .groupby("id")
    selected = []
    for name, group in affgrp:
        wt_aff = group[group["comment"] == "wt"].iloc[0]["%s_score" % sitetarget]
        wt_pred = group[group["comment"] == "wt"].iloc[0]["main_pred"]
        mt_selected = group[(group["comment"] == "wt") |
                    ((group["main_pred"] != wt_pred) & (abs(wt_aff - group["%s_score" % sitetarget]) < affdelta))]
        if mt_selected.shape[0] > 1:
            selected.extend(mt_selected.to_dict('records'))
    select_df = pd.DataFrame(selected)
    select_df["select"] = "small_affdif_%s" % sitetarget
    return select_df

def count_selected(df):
    count = df["muttype"].value_counts()
    print("Count muttype", count)

np.random.seed(100)
pd.set_option("display.max_columns",None)
if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"

    filename = "output/mutall_etsrunx_delta.csv"
    df = pd.read_csv(filename)
    df = df[~df['Sequence'].str.contains('GGGGG') & ~df['Sequence'].str.contains('CCCCC')]
    get_count(df)

    # get entries where wt predictions match
    matchwt = df[(df["comment"] == "wt") & (df["wtlabel"] == df["main_pred"])].drop_duplicates("Sequence")
    matchwt_ids = matchwt[["id"]].drop_duplicates()
    df = df.merge(matchwt_ids, on="id").sort_values(by=["id","muttype","distance","s1_score","s2_score"], ascending=[True,True,False,True,True])

    muttypes = ["strength", "distance", "orientation"]

    # initialize with wt
    allpicks = pd.DataFrame(matchwt)
    allpicks["select"] = "wt"

    # allpicks.to_csv("bla.csv",index=False)

    for mty in muttypes:
        curdf = df.loc[df["muttype"] == mty]
        posctrl = pick_positive_ctrl(curdf, 200, 40, "main", "shape")
        p1 = pick_largest_probdif(curdf, 200, 40, mty)
        p2 = pick_incos_pred(curdf, 140, 40)
        p = pd.concat([p1, p2, posctrl])
        if mty == "strength":
            p_s1 = pick_small_affdif(curdf, 120, 30, "s1")
            p_s2 = pick_small_affdif(curdf, 120, 30, "s2")
            p = pd.concat([p,p_s1,p_s2])
        p = p[p["comment"] != "wt"]
        allpicks = pd.concat([allpicks, p])

    allpicks = allpicks[~allpicks['Sequence'].str.contains('GGGGG') & ~allpicks['Sequence'].str.contains('CCCCC')]
    allpicks = allpicks[allpicks['id'].groupby(allpicks['id']).transform('size')>1]
    suffix = "etsrunx" if "etsrunx" in filename else "etsets"
    allpicks.sort_values(["id","comment"]).to_csv("custom_selected_%s.csv" % suffix, index=False, header=True)

    gap = allpicks.groupby("id")
    for a,b in gap:
        x = b[b["comment"] == "wt"]
        if x.empty:
            print("empty",a)

    count_selected(allpicks)
    print("Number of unique sequences %d" % len(allpicks["Sequence"].unique()))
    totalspots = 16 if "etsrunx" in filename else 24
    print("Total array spots",len(allpicks["Sequence"].unique()) * totalspots)
