import pandas as pd
import math
import itertools
import pickle

from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.sitesplotter import SitesPlotter
import chip2probe.training_gen.traingen as tg
from chip2probe.util import bio

def eliminate_site(seq, kompdict, pwmdict, target, maxpwmdelta=1):
    bases = {'A', 'T', 'G', 'C'}
    pred = tg.pred_all(seq, kompdict, pwmdict)
    orig_sites_count = len(pred)
    pred_core_starts = [p["core_start"] for p in pred]
    lowest_score = pred[target]["score"] # assign lowest score with the target score
    targeted_core_start = pred[target]["core_start"]
    untargeted_core_start = {pred[i]["core_start"]:pred[i]["score"] for i in range(len(pred)) if i != target}

    pred = pred[target]
    pwmtf = pwmdict[pred['tf']]
    flanklen = (pwmtf.length - pred["core_width"])//2
    s_start, s_end = pred["core_start"]-flanklen,pred["core_start"]+pred["core_width"]+flanklen

    lbase, lpos = "", -1
    mutpos = list(range(pred["core_start"], pred["core_start"]+pred["core_width"]))
    saved_mutpos = []
    site_exist = True
    mutating = True
    mutseq = str(seq)
    while site_exist and mutating:
        for pos in mutpos:
            for b in bases:
                curmut = seq[:pos] + b + seq[pos+1:]
                mutscr = pwmtf.predict_sequence(curmut[s_start:s_end], zero_thres=False)[0]["score"]
                if mutscr < lowest_score:
                    curmutseq = seq[:pos] + b + seq[pos+1:]
                    # make sure the mutation doesn't create new site nor change the strength of the other site
                    curmutpred = tg.pred_all(curmutseq, kompdict, pwmdict)
                    mut_ok = True
                    if len(curmutpred) < orig_sites_count-1:
                        # we destroyed the other site
                        mut_ok = False
                    else:
                        for i in range(len(curmutpred)):
                            if curmutpred[i]["core_start"] == targeted_core_start:
                                # current site
                                continue
                            if not curmutpred[i]["core_start"] in pred_core_starts:
                                # we created new sites
                                mut_ok = False
                                break
                            if not math.isclose(curmutpred[i]["score"],untargeted_core_start[curmutpred[i]["core_start"]],rel_tol=maxpwmdelta):
                                # we changed an existing site
                                mut_ok = False
                                break
                    if mut_ok:
                        lowest_score = mutscr
                        lbase, lpos = b, pos
        mutseq = seq[:lpos] + lbase + seq[lpos+1:]
        mutpred = tg.pred_all(mutseq, kompdict, pwmdict)
        if not targeted_core_start in [m["core_start"] for m in mutpred]:
            site_exist = False
        # we have targeted this position before
        if lpos not in saved_mutpos:
            saved_mutpos.append(lpos)
        else:
            # we have made the best mutation
            mutating = False
    if not mutating:
        return ""
    else:
        return mutseq

def make_mutations(seqlist, kompdict, pwmdict, toelims=[0,1,2]):
    mutlist = []
    md = len(seqlist) // 100 if len(seqlist) > 100 else -1
    for i in range(len(seqlist)):
        if md != -1 and i % md == 0:
            print("Progress mutating %d/%d" % (i+1,len(seqlist)))
        seq = seqlist[i]
        if len(tg.pred_all(seq, kompdict, pwmdict)) != 2:
            print("Mutating error for sequence %s" % seq)
            continue
        mutating = True
        mutdict = {x:"" for x in toelims}
        for target in toelims:
            if target != 2:
                mut = eliminate_site(seq, kompdict, pwmdict, target)
            else:
                mut = eliminate_site(mutdict[0], kompdict, pwmdict, 0)
            if mut != "":
                mutdict[target] = mut
            else:
                mutating = False
                break
        if mutating:
            mutd = {"wt":seq, "m1":mutdict[0] ,"m2":mutdict[1]}
            if len(toelims) > 2:
                mutd["m3"] = mutdict[2]
            mutlist.append(mutd)
        else:
            print("can't mutate %s" % seq)
    return mutlist

def is_junction_site(seq, primer, kompdict, pwmdict):
    rc = bio.revcompstr(seq)
    sqp = seq + primer
    rcp = rc + primer
    sqp_pred = tg.pred_all(sqp, kompdict, pwmdict)
    rcp_pred = tg.pred_all(rcp, kompdict, pwmdict)
    if len(sqp_pred) > 2:
        return "fwd"
    elif len(rcp_pred) > 2:
        return "rev"
    else:
        return "ok"

def clean_junction(seqlist, primer, kompdict, pwmdict):
    allseqs = []
    modifiedseqs = []
    md = len(seqlist) // 100 if len(seqlist) > 100 else -1
    for i in range(len(seqlist)):
        if md != -1 and i % md == 0:
            print("Progress %d/%d" % (i+1,len(seqlist)))
        seq = seqlist[i]
        seqlen = len(seq)
        stats = is_junction_site(seq, primer, kompdict, pwmdict)
        if stats != "ok":
            targetseq = str(seq) if stats == "fwd" else bio.revcompstr(seq)
            newseq = eliminate_site(targetseq+primer, kompdict, pwmdict, 2)[:seqlen]
            if is_junction_site(newseq, primer, kompdict, pwmdict) != "ok":
                print("can't mutate")
                continue
            else:
                newseq = bio.revcompstr(newseq) if stats == "rev" else newseq
                modifiedseqs.append({"old": seq,"new":newseq})
                allseqs.append(newseq)
        else:
            allseqs.append(seq)
    return allseqs, modifiedseqs

def make_custom(seqlist, primer, prefix, rep=3, arrname="Coop3Ets_"):
    arrstr = ""
    id = 0
    maxidlen = 15
    zf = maxidlen - len(arrname)

    for i in range(len(seqlist)):
        for o in ["o1", "o2"]:
            for key in seqlist[i]:
                if o == "o2":
                    curseq = bio.revcompstr(seqlist[i][key]) + primer
                else:
                    curseq = seqlist[i][key] + primer
                for j in range(1,rep+1):
                    id += 1
                    cur_id = "%s%s" % (arrname, str(id).zfill(zf))
                    type = "%s_%d_%s_%s_r%d" % (prefix, i+1, key, o, j)
                    arrstr += "%s\t%s\t%s\tNa|Na\tNa\t%s\tchr1:0-0\n" % (cur_id, curseq, type, type)
    return arrstr

pd.set_option("display.max_columns",None)
if __name__ == "__main__":
    primer = "GTCTTGATTCGCTTGACGCTGCTG"
    filepath = "output/custom_selected_etsrunx.csv"
    # filepath = "input/litseqs_etsets.csv"

    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"
    pwm_ets = PWM("%s/input/sitemodels/pwm/ets1.txt" % basepath, log=True, reverse=False)
    kompas_ets = Kompas("%s/input/sitemodels/kompas/Ets1_kmer_alignment.txt" % basepath, core_start = 11, core_end = 15, core_center = 12)
    pwm_runx = PWM("%s/input/sitemodels/pwm/runx1.txt" % basepath, 8, 17, log=True, reverse=True)
    kompas_runx = Kompas("%s/input/sitemodels/kompas/Runx1_kmer_alignment.txt" % basepath, core_start = 12, core_end = 17, core_center = 14)

    kompdict = {"ets1":kompas_ets}
    pwmdict = {"ets1":pwm_ets}
    if "etsrunx" in filepath:
        toelims = [0,1]
        kompdict["runx1"] = kompas_runx
        pwmdict["runx1"] = pwm_runx
        prefix = "er"
        filesuf = "etsrunx"
    else:
        toelims = [0,1,2]
        prefix = "ee"
        filesuf = "etsets"

    df = pd.read_csv(filepath)
    allseqs = df["Sequence"].unique().tolist()
    print("Working on %d sequences" % len(allseqs))

    allseqs_junc, modseqs = clean_junction(allseqs, primer, kompdict, pwmdict)
    pd.DataFrame(modseqs).to_csv("cleanjunction_map_%s.csv" % filesuf, index=False)
    muts = make_mutations(allseqs_junc, kompdict, pwmdict, toelims=toelims)

    with open("arr_%s.txt" % filesuf,"w") as f:
        f.write(make_custom(muts, primer, prefix))

    # allseqs = [elm[k] for elm in muts for k in elm]
    # kompas_plot = [kompdict[k].make_plot_data(kompdict[k].predict_sequences(allseqs)) for k in kompdict]
    # sp = SitesPlotter()
    # sp.plot_seq_combine(kompas_plot, filepath="plot_%s.pdf"%filesuf,numcol=len(toelims)+1)
