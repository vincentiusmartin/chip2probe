'''
Created on May, 2020

@author: vincentiusmartin
@description: Make mutations for custom probes. All columns on the mutants are
    repredicted so we can validate that mutations don't change features that
    should not be changed.
'''

import chip2probe.util.coopgeneral as cg
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.util import bio

import math
import pandas as pd
import itertools

def mutate_affinity(seqdf, imads, deep=0):
    """
    Make mutation to change the affinity (i.e. strength) prediction.

    First, mutations are made for each core to its other core versions, e.g. if
    the core is GGAA and the alternate is GGAT, we simply change GGAAA -> GGAT.
    Then mutate the core flanking regions up to imads.sitewidth, e.g. if the
    core length is 4 and sitewidth is 12 then we can mutate up to (12-4)/2=4bp
    to each side. When 'deep' is set to be more than 0, set barrier to
    sitewidth - distance on the other binding site.

    Args:
        1. seqdf: input data frame with the wt sequences to mutate
        2. imads: imads model to predict the strength of the mutants
        3. deep: the minimum distance between sequence is set to be
            imads.sitewidth - deep. Default is 0, which means we keep sitewidth
            as minimum distance. When deep is > 0, we make barrier at the other
            site so we don't change its affinity prediction.
     Returns:
        A data frame of sequences with SNPs that change its affinity.
    """
    if deep < 0:
        raise ValueError("Minimum deep is 0")

    ct = CoopTrain(seqdf["sequence"].values.tolist(), corelen=4, flip_th=True, imads=imads, ignore_sites_err=True)
    om = ct.df.join(seqdf.set_index("sequence"), on="sequence", how="inner") # this already include the orientation

    # first make map for mutating between core
    mdlcores_fw = [m.core for m in imads.models]
    fwdict = {e[0]:e[1] for e in list(itertools.permutations(mdlcores_fw,2))}
    mdlcores_rc = [bio.revcompstr(m) for m in mdlcores_fw]
    rcdict = {e[0]:e[1] for e in list(itertools.permutations(mdlcores_rc,2))}
    coremap = {**fwdict, **rcdict}

    # prepare the variable
    mindist = imads.sitewidth - deep
    mutres = []

    iter = 0
    nrow = om.shape[0]
    div = nrow // 100
    for index, row in om.iterrows():
        if iter % div == 0:
            print("Mutating affinity, progress {:.2f}% ({}/{})".format(iter*100/nrow,iter,nrow))
        iter += 1
        mutres_cur = []
        sites = imads.predict_sequence(row["sequence"])
        if sites[1]["core_mid"] - sites[0]["core_mid"] < mindist or len(sites) != 2:
            continue
        mutres_cur.append({
            "seqid":row["id"],
            "sequence":str(row["sequence"]),
            "site1_pos":sites[0]["core_mid"],
            "site1_affinity":sites[0]["score"],
            "site2_pos":sites[1]["core_mid"],
            "site2_affinity":sites[1]["score"],
            "distance":sites[1]["core_mid"] - sites[0]["core_mid"],
            "muttype":"affinity",
            "comment":"wt",
            "wtlabel":row["label"],
            "orientation":row["orientation"]
        })
        mids = [s["core_mid"] for s in sites]

        # 1. Mutate the core to the other version
        coremt = mutate_cores(row["sequence"], mids, coremap)

        # 2. Mutate the flanks up to the sitewidth
        barrierlen = imads.sitewidth - row["distance"] if row["distance"] < imads.sitewidth else 0
        flankmt = mutate_flanks(row["sequence"], mids, imads.corewidth, imads.sitewidth, barrier=barrierlen)

        allmt = coremt + flankmt
        for i in range(len(allmt)):
            newsites = imads.predict_sequence(allmt[i]["sequence"])
            if len(newsites) != 2:
                continue
            newori = cg.get_relative_orientation(allmt[i]["sequence"], imads, htth=True)
            mutres_cur.append({
                "seqid":row["id"],
                "sequence":allmt[i]["sequence"],
                "site1_pos":newsites[0]["core_mid"],
                "site1_affinity":newsites[0]["score"],
                "site2_pos":newsites[1]["core_mid"],
                "site2_affinity":newsites[1]["score"],
                "distance":newsites[1]["core_mid"]-newsites[0]["core_mid"],
                "muttype":"affinity",
                "comment":allmt[i]["comment"],
                "wtlabel":row["label"],
                "orientation":newori
            })
        if len(mutres_cur) > 1:
            mutres.extend(mutres_cur)
    return pd.DataFrame(mutres)

def mutate_flanks(seq, midlist, corelen, sitelen, barrier=0):
    atgc = {'A', 'T', 'G', 'C'}
    coreside = corelen // 2
    flank = sitelen // 2 - coreside
    mutants = []
    for i in range(len(midlist)):
        mid = midlist[i]
        lpos = list(range(mid-coreside-1, mid-coreside-flank-1,-1))
        rpos = list(range(mid+coreside, mid+coreside+flank))
        poslist = lpos + rpos
        posskip = list(range(mid+coreside+flank-1, mid+coreside+flank-barrier-1, -1)) if i == 0 \
                  else list(range(mid-coreside-flank, mid-coreside-flank + barrier))
        for p in poslist:
            if p in posskip:
                continue
            lr = "left" if p < mid else "right"
            lrpos = mid-p-2 if p < mid else p-mid-1
            for nuc in atgc - set(seq[p]):
                comment = "site%d_%s_pos%d_%s" % (i+1, lr, lrpos, nuc)
                mutants.append({"sequence":seq[:p] + nuc + seq[p + 1:],
                                "comment": comment})
    return mutants

def mutate_cores(seq, midlist, coremap):
    # assume all cores have the same length
    corelen = len(list(coremap.keys())[0])
    sidelen = corelen // 2
    mutants = []
    for i in range(len(midlist)):
        mid = midlist[i]
        start, end = mid-sidelen, mid+sidelen
        core = seq[start:end]
        comment = "core_site%s" % (i+1)
        mutants.append({"sequence":seq[:start] + coremap[core] + seq[end:],
                        "comment": comment})
    return mutants

def mutate_orientation(seqdf, imads, deep=0):
    """
    Make mutation for orientation

    Flip one or both sites, flip the whole 12 mer. We only use HH, HT, TT
    orientation (i.e. no TH).

    Args:
        1. seqdf: input data frame with the wt sequences to mutate
        2. imads: imads model to predict the strength of the mutants
        3. deep: how far we permit distance to go under imads.sitewidth. The
            minimum distance is set to imads.sitewidth - deep. The flip length
            is changed from sitewidth to (sitewidth-deep)//2*2.
     Returns:
        A data frame with changed orientations
    """
    # we need to get orientation information, this already filter if each sequence has 2 sites
    ct = CoopTrain(seqdf["sequence"].values.tolist(), corelen=4, flip_th=True, imads=imads, ignore_sites_err=True)
    om = ct.df.join(seqdf.set_index("sequence"), on="sequence", how="inner") # this already include the orientation
    mutres = []
    orilist = {"HH","TT","HT/TH"}
    flipsites = [[0],[1],[0,1]] # which sites to flip
    iter = 0
    nrow = om.shape[0]
    div = nrow // 100
    mindist = imads.sitewidth - deep
    for index, row in om.iterrows():
        if iter % div == 0:
            print("Mutating orientation, progress {:.2f}% ({}/{})".format(iter*100/nrow,iter,nrow))
        iter += 1
        mutres_cur = []
        sites = imads.predict_sequence(row["sequence"])
        curdist = sites[1]["core_mid"] - sites[0]["core_mid"]
        if  curdist < mindist or len(sites) != 2:
            continue
        mutres_cur.append({
            "seqid":row["id"],
            "sequence":str(row["sequence"]),
            "site1_pos":sites[0]["core_mid"],
            "site1_affinity":sites[0]["score"],
            "site2_pos":sites[1]["core_mid"],
            "site2_affinity":sites[1]["score"],
            "distance":sites[1]["core_mid"]-sites[0]["core_mid"],
            "muttype":"orientation",
            "comment":"wt",
            "wtlabel":row["label"],
            "orientation":row["orientation"]
        })
        for fs in flipsites:
            newseq = row["sequence"]
            adjust = 0 if curdist >= imads.sitewidth else int(math.ceil(float(imads.sitewidth - curdist) / 2))
            for i in fs:
                start,end = sites[i]["site_start"]+adjust, sites[i]["site_start"]+sites[i]["site_width"]-adjust
                toflip = bio.revcompstr(row["sequence"][start:end])
                newseq = newseq[:start] + toflip + newseq[end:]
            newsites = imads.predict_sequence(newseq)
            if len(newsites) != 2: # we ignore if there are new sites
                continue
            newori = cg.get_relative_orientation(newseq, imads, htth=False)
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
        # if len(mutres_cur) != 3: # 3 orientations
        #     print("Found predictions with number of orientation probes != 3", len(mutres_cur), row["sequence"])
        if len(mutres_cur) > 1:
            mutres.extend(mutres_cur)
    return pd.DataFrame(mutres)

def mutate_dist(seqdf, imads, deep=0,  warning=True):
    """
    Make mutation for distance

    Make closer distance between two sites.
    Insert and cut, for cutting can fix the second site (site that closer to the
    glass slide) and can just use the the nucleotide that is cut to patch.
    CHECK WE ARE NOT CREATING NEW SITE.

    Args:
        feature_dict: the following is the list of currently available feature:
            1. seqdf: input data frame with the wt sequences to mutate
            2. imads: imads model to predict the strength of the mutants
            3. deep: how far we permit distance to go under imads.sitewidth. The
                minimum distance is set to imads.sitewidth - deep. When deep > 0,
                the affinity of each site will change after the distance is
                less than sitewidth.
            4. warning: print warning when input has sites number != 2
    Returns:
        A data frame where each sequence has mutants
    """
    ct = CoopTrain(seqdf["sequence"].values.tolist(), corelen=4, flip_th=True, imads=imads, ignore_sites_err=True)
    om = ct.df.join(seqdf.set_index("sequence"), on="sequence", how="inner") # this already include the orientation

    mutres = []
    mindist = imads.sitewidth - deep
    nrow = om.shape[0]
    iter = 0
    div = nrow // 100
    for index, row in om.iterrows():
        if iter % div == 0:
            print("Mutating distance, progress {:.2f}% ({}/{})".format(iter*100/nrow,iter,nrow))
        iter += 1
        seq = row["sequence"]
        mutres_cur = []
        sites = imads.predict_sequence(seq)
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
            "muttype":"distance",
            "comment":"wt",
            "wtlabel":row["label"],
            "orientation":row["orientation"]
        })
        move = 1
        initdist = int(curdist) # need to make this because we keep changing curdist
        while initdist-move >= mindist:
            s1_end = sites[0]["site_start"]+sites[0]["site_width"]
            s2_start = sites[1]["site_start"]
            curseq = cg.move_single_site(seq, s1_end, s2_start, move, patch=True, can_overlap=True)
            cursites = imads.predict_sequence(curseq)
            move += 1
            if len(cursites) != 2: # we ignore if there are new sites
                continue
            curdist = cursites[1]["core_mid"] - cursites[0]["core_mid"]
            newori = cg.get_relative_orientation(curseq, imads, htth=True)
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
                "wtlabel":row["label"],
                "orientation":newori
            })
        if len(mutres_cur) > 1:
            mutres.extend(mutres_cur)
    return pd.DataFrame(mutres)
