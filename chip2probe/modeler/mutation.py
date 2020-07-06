import pandas as pd
import chip2probe.util.coopgeneral as cg

# 1446
def pick_custom_dist(distdf):
    grouped = distdf.groupby('seqid')
    shiftedpred = []
    shiftcount = 0
    startcoopct = 0
    startaddct = 0
    wrongpred_count = 0
    #print(len(grouped))
    for name, group in grouped:
        first = group.head(1).iloc[0]
        # first, ignore if prediction is different with wt
        if first["wtlabel"] != first["main_pred"]:
            wrongpred_count += 1
            continue
        last = group.tail(1).iloc[0]
        if first["main_pred"] != last["main_pred"]:
            print(first["main_pred"])
            if first["main_pred"] == 1:
                startcoopct += 1
            else:
                startaddct += 1
            shiftedpred.extend(group.to_dict('records'))
            shiftcount += 1
    print(startaddct, startcoopct)
    pd.DataFrame(shiftedpred).to_csv("sp.csv", index=False, header=True)

def mutate_dist(seqdf, predictor, fixsecond=True, warning=True):
    """
    Make mutation for distance

    Insert and cut, for cutting can fix the second site (site that closer to the glass slide) and can just use the the nucleotide that is cut to patch. CHECK WE ARE NOT CREATING NEW SITE.

    Args:
        feature_dict: the following is the list of currently available feature:
            1. df:
            2: pred:
            3: fixsecond
     Returns:

    """
    mutres = []
    mindist = predictor.sitewidth
    for index, row in seqdf.iterrows():
        seq = row["sequence"]
        mutres_cur = []
        sites = predictor.predict_sequence(seq)
        if len(sites) != 2:
            if warning:
                print("Found a sequence with number of sites not equal to 2: ",seq,sites)
            continue
        curdist = sites[1]["core_mid"] - sites[0]["core_mid"]
        if curdist <= mindist: # if less or equal, we can't mutate
            continue
        mutres_cur.append({
            "sequence":str(seq),
            "seqid":"seq%d" % index,
            "site1_pos":sites[0]["core_mid"],
            "site1_affinity":sites[0]["score"],
            "site2_pos":sites[1]["core_mid"],
            "site2_affinity":sites[1]["score"],
            "distance":int(curdist),
            "comment":"wt",
            "wtlabel":row["label"]
        })
        move = 1
        while curdist > mindist:
            s1_end = sites[0]["site_start"]+sites[0]["site_width"]
            s2_start = sites[1]["site_start"]
            curseq = cg.move_single_site(seq, s1_end, s2_start, move, patch=True)
            cursites = predictor.predict_sequence(curseq)
            if len(cursites) != 2: # we stop iterating
                break
            curdist = cursites[1]["core_mid"] - cursites[0]["core_mid"]
            mutres_cur.append({
                "sequence":str(curseq),
                "seqid":"seq%d" % index,
                "site1_pos":cursites[0]["core_mid"],
                "site1_affinity":cursites[0]["score"],
                "site2_pos":cursites[1]["core_mid"],
                "site2_affinity":cursites[1]["score"],
                "distance":int(curdist),
                "comment":"closer_%d" % move,
                "wtlabel":row["label"]
            })
            move += 1
        if len(mutres_cur) > 1:
            mutres.extend(mutres_cur)
    return pd.DataFrame(mutres)
