import pandas as pd



def mutate_dist(sequences, predictor, fixsecond=True, warning=True):
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
    for seq in sequences:
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
            "site_str_pos":sites[0]["core_mid"] if sites[0]["score"] > sites[1]["score"] else sites[1]["core_mid"],
            "site_str_affinity":sites[0]["score"] if sites[0]["score"] > sites[1]["score"] else sites[1]["score"],
            "site_wk_pos":sites[0]["core_mid"] if sites[0]["score"] < sites[1]["score"] else sites[1]["core_mid"],
            "site_wk_affinity":sites[0]["score"] if sites[0]["score"] < sites[1]["score"] else sites[1]["score"],
            "distance":int(curdist),
            "comment":"wt"
        })
        move = 1
        while curdist > mindist:
            s1b = sites[0]["site_start"]+sites[0]["site_width"]
            s2a = sites[1]["site_start"]
            replaced = seq[s1b:s1b+move]
            curseq = replaced + seq[:s1b] + seq[s1b+move:]
            cursites = predictor.predict_sequence(curseq)
            if len(cursites) != 2: # we stop iterating
                break
            curdist = cursites[1]["core_mid"] - cursites[0]["core_mid"]
            mutres_cur.append({
                "sequence":str(curseq),
                "site_str_pos":cursites[0]["core_mid"] if cursites[0]["score"] > cursites[1]["score"] else cursites[1]["core_mid"],
                "site_str_affinity":cursites[0]["score"] if cursites[0]["score"] > cursites[1]["score"] else cursites[1]["score"],
                "site_wk_pos":cursites[0]["core_mid"] if cursites[0]["score"] < cursites[1]["score"] else cursites[1]["core_mid"],
                "site_wk_affinity":cursites[0]["score"] if cursites[0]["score"] < cursites[1]["score"] else cursites[1]["score"],
                "distance":int(curdist),
                "comment":"closer_%d" % move
            })
            move += 1
        if len(mutres_cur) > 1:
            mutres.extend(mutres_cur)
    return pd.DataFrame(mutres)
