import pandas as pd
import chip2probe.util.coopgeneral as cg
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.util import bio

import pandas as pd

def mutate_orientation(seqdf, predictor):
    """
    Make mutation for orientation

    Flip one or both sites, flip the whole 12 mer

    Args:
        feature_dict: the following is the list of currently available feature:
            1. seqdf:
            2: predictor
     Returns:

    """
    # we need to get orientation information, this already filter if each sequence has 2 sites
    ct = CoopTrain(seqdf["sequence"].values.tolist(), corelen=4, flip_th=True, imads=predictor, ignore_sites_err=True)
    om = ct.df.join(seqdf.set_index("sequence"), on="sequence", how="inner") # this already include the orientation
    mutres = []
    orilist = {"HH","TT","HT/TH"}
    flipsites = [[0],[1],[0,1]] # which sites to flip
    iter = 0
    nrow = om.shape[0]
    div = nrow // 100
    for index, row in om.iterrows():
        if iter % div == 0:
            print("Mutating orientation, progress {:.2f}% ({}/{})".format(iter*100/nrow,iter,nrow))
        iter += 1
        mutres_cur = []
        sites = predictor.predict_sequence(row["sequence"])
        mutres_cur.append({
            "seqid":row["id"],
            "sequence":str(row["sequence"]),
            "site1_pos":sites[0]["core_mid"],
            "site1_affinity":sites[0]["score"],
            "site2_pos":sites[1]["core_mid"],
            "site2_affinity":sites[1]["score"],
            "distance":row["distance"],
            "muttype":"wt",
            "comment":"wt",
            "wtlabel":row["label"],
            "orientation":row["orientation"]
        })
        sites = predictor.predict_sequence(row["sequence"])
        for fs in flipsites:
            newseq = row["sequence"]
            for i in fs:
                start,end = sites[i]["site_start"], sites[i]["site_start"]+sites[i]["site_width"]
                toflip = bio.revcompstr(row["sequence"][start:end])
                newseq = newseq[:start] + toflip + newseq[end:]
            newsites = predictor.predict_sequence(newseq)
            if len(newsites) != 2: # we ignore if there are new sites
                continue
            newori = cg.get_relative_orientation(newseq, predictor, htth=False)
            if newori == "HT":
                newori = "HT/TH"
            elif newori == "TH":
                continue # skip if TH since we use HT
            mutres_cur.append({
                "seqid":row["id"],
                "sequence":str(newseq),
                "site1_pos":newsites[0]["core_mid"],
                "site1_affinity":newsites[0]["score"],
                "site2_pos":newsites[1]["core_mid"],
                "site2_affinity":newsites[1]["score"],
                "distance":newsites[1]["core_mid"]-newsites[0]["core_mid"],
                "muttype":"orientation",
                "comment":"to_%s" % newori,
                "wtlabel":row["label"],
                "orientation":newori
            })
        if len(mutres_cur) > 1:
            mutres.extend(mutres_cur)
    return pd.DataFrame(mutres)

def mutate_dist(seqdf, predictor, fixsecond=True, warning=True):
    """
    Make mutation for distance

    Insert and cut, for cutting can fix the second site (site that closer to the glass slide) and can just use the the nucleotide that is cut to patch. CHECK WE ARE NOT CREATING NEW SITE.

    Args:
        feature_dict: the following is the list of currently available feature:
            1. df: input data frame with id, sequence, label
            2: predictor:
            3: fixsecond
     Returns:

    """
    mutres = []
    mindist = predictor.sitewidth
    nrow = seqdf.shape[0]
    iter = 0
    div = nrow // 100
    for index, row in seqdf.iterrows():
        if iter % div == 0:
            print("Mutating distance, progress {:.2f}% ({}/{})".format(iter*100/nrow,iter,nrow))
        iter += 1
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
            "seqid":row["id"],
            "sequence":str(seq),
            "site1_pos":sites[0]["core_mid"],
            "site1_affinity":sites[0]["score"],
            "site2_pos":sites[1]["core_mid"],
            "site2_affinity":sites[1]["score"],
            "distance":int(curdist),
            "muttype":"wt",
            "comment":"wt",
            "wtlabel":row["label"]
        })
        move = 1
        initdist = int(curdist)
        while initdist-move >= mindist:
            s1_end = sites[0]["site_start"]+sites[0]["site_width"]
            s2_start = sites[1]["site_start"]
            curseq = cg.move_single_site(seq, s1_end, s2_start, move, patch=True)
            cursites = predictor.predict_sequence(curseq)
            move += 1
            if len(cursites) != 2: # we ignore if there are new sites
                continue
            curdist = cursites[1]["core_mid"] - cursites[0]["core_mid"]
            mutres_cur.append({
                "seqid":row["id"],
                "sequence":str(curseq),
                "site1_pos":cursites[0]["core_mid"],
                "site1_affinity":cursites[0]["score"],
                "site2_pos":cursites[1]["core_mid"],
                "site2_affinity":cursites[1]["score"],
                "distance":int(curdist),
                "muttype":"distance",
                "comment":"closer_%d" % (move-1),
                "wtlabel":row["label"]
            })
        if len(mutres_cur) > 1:
            mutres.extend(mutres_cur)
    return pd.DataFrame(mutres)
