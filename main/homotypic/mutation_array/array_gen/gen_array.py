import os
os.chdir("../../../..")
import pandas as pd
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pbmescore import PBMEscore
from chip2probe.sitespredict.dnasequence import DNASequence
import chip2probe.util.bio as bio

def sanity_check(df, imads, primer):
    print("sanity check...")
    invalid_ids = []
    i = 0
    for index, row in df.iterrows():
        if i % 300 == 0:
            print("Progress %d/%d" % (i+1,df.shape[0]))
        i += 1
        nsites = len(imads.predict_sequence(row["sequence"]))
        # if not len(row["sequence"]) == 60:
        #     raise Exception("Wrong sequence length", row["id"])
        if not row["sequence"].endswith(primer):
            raise Exception("primer is missing", row["id"])
        if ("wt" in row["id"] and nsites != 2) or \
           (("m1" in row["id"] or "m2" in row["id"])  and nsites != 1) or \
           ("m3" in row["id"]  and nsites != 0) or \
           ("NegativeCtrl" in row["id"] and nsites != 0):
            raise Exception("wrong site count", row["sequence"], row["id"],nsites)
        #invalid_ids.append()
    return True
# check primer exist

def clean_junction(sequence, primer, imads, escore):
    s_site = imads.predict_sequence(sequence)
    sp = sequence+primer
    sp_site = imads.predict_sequence(sp)
    if len(s_site) == len(sp_site):
        # no new site, we are good
        return sp
    else:
        ds = DNASequence(sp, imads12, escore)
        to_mutate = len(s_site)
        newseq = ds.abolish_sites([to_mutate],nomutrange=[36,len(sp)])
        return newseq

if __name__ == "__main__":
    primer = "GTCTTGATTCGCTTGACGCTGCTG"
    arrayname = "Coop2Ets_"
    maxidlen = 15
    maxprobenum = 176024
    zf = maxidlen - len(arrayname) # zero fill

    # Load escore object
    escore = PBMEscore("input/site_models/escores/Ets1_8mers_11111111.txt")

    # Load imads object
    imads12_paths = ["input/site_models/imads_models/Ets1_w12_GGAA.model", "input/site_models/imads_models/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19) # 0.2128

    df = pd.read_csv("arrayseqs.csv")
    indict = pd.Series(df["sequence"].values,index=df["id"]).to_dict()
    seqids = df.loc[df["id"].str.contains("wt")]["id"].tolist()
    seqids = [x.split("_")[0] for x in seqids]
    arr = []

    negdf = pd.read_csv("output/array_design_files/Coop1Ets/Coop1Ets_NegCtrl.txt",sep="\t",header=None)
    print(negdf)

    negctrl = negdf[0].tolist()
    negids = [x.split("_")[3] for x in negdf[2].tolist()]
    validneg = []
    for neg_idx in range(0,len(negctrl)):
        cur_pos_seq = negctrl[neg_idx]
        for ori in ["o1", "o2"]:
            checkseq = cur_pos_seq if ori == "o1" else bio.revcompstr(cur_pos_seq)
            cj = clean_junction(checkseq, primer, imads12, escore)
            if not cj:
                print("Problem in clean junction")
                flag = False
                break
            if cj != checkseq+primer:
                if ori == "o1":
                    cur_pos_seq = str(cj[:len(cur_pos_seq)])
                elif ori == "o2": # change o1
                    cur_pos_seq = bio.revcompstr(str(cj[:len(cur_pos_seq)]))
        if len(imads12.predict_sequence(cur_pos_seq+primer)) > 0:
            print("Error processing sequence",cur_pos_seq)
        else:
            validneg.append({"sequence":cur_pos_seq,"type":negids[neg_idx]})

    maxwtnum = (maxprobenum - (len(validneg) * 6)) // 24

    vtype = ["wt","m1","m2","m3"]
    wtlist = []
    for i in range(len(seqids)):
        if i % 70 == 0:
            print("Progress %d/%d" % (i+1,len(seqids)))
        cur_arr = []
        flag = True
        # do clean junction first on the wt_o1
        o1_wt = indict["%s_wt"%(seqids[i])]
        cleanedseq = indict["%s_wt"%(seqids[i])]
        for ori in ["o1", "o2"]:
            curseq = o1_wt if ori == "o1" else bio.revcompstr(o1_wt)
            cj = clean_junction(curseq, primer, imads12, escore)
            if not cj:
                print("Problem in clean junction")
                flag = False
                break
            if cj != curseq+primer:
                if ori == "o1":
                    cleanedseq = str(cj[:len(cleanedseq)])
                elif ori == "o2": # change o1
                    cleanedseq = bio.revcompstr(str(cj[:len(cleanedseq)]))
        if cleanedseq != indict["%s_wt"%(seqids[i])]:
            # we change the sequence so we need to update them all
            posdif = [i for i in range(len(curseq)) if o1_wt[i] != cleanedseq[i]]
            for v in vtype:
                for p in posdif:
                    orig = indict["%s_%s"%(seqids[i],v)]
                    indict["%s_%s"%(seqids[i],v)] = orig[:p] + cleanedseq[p] + orig[p+1:]

        for ori in ["o1","o2"]:
            for v in vtype:
                s = indict["%s_%s"%(seqids[i],v)]
                wprimer = s if ori == "o1" else bio.revcompstr(s)
                wprimer += primer
                nsites = len(imads12.predict_sequence(wprimer))
                if ("wt" in v and nsites != 2) or \
                   (("m1" in v or "m2" in v)  and nsites != 1) or \
                   ("m3" in v  and nsites != 0):
                    print("wrong site count", wprimer, v,nsites)
                    flag = False
                    break
                if not wprimer.endswith(primer): # or len(wprimer) != 60:
                    raise Exception("Problematic sequence")
                cur_arr.append({
                    "sequence": wprimer,
                    "id": "%s_%s" % (seqids[i],v),
                    "ori": ori
                })
            if not flag:
                print("Error processing sequence %d" % i)
                break # ignore
        if flag and indict["%s_wt"%(seqids[i])] not in wtlist:
            wtlist.append(indict["%s_wt"%(seqids[i])])
            arr.extend(cur_arr)
            if len(wtlist) >= maxwtnum:
                break
    pd.DataFrame(arr).to_csv("finarr.csv",index=False, header=True)
    # # cur_id = "%s%s" % (arrayname, str(i+1).zfill(zf))
    df = pd.read_csv("finarr.csv")

    rowlist = []
    i = 1
    for index,row in df.iterrows():
        for r in ["r1","r2","r3"]:
            cur_id = "%s%s" % (arrayname, str(i).zfill(zf))
            rowdesc = "%s_%s_%s" % (row["id"],row["ori"],r)
            currow = "%s\t%s\t%s\tNa|Na\tNa\t%s\tchr1:0-0" % (cur_id, row["sequence"], rowdesc , rowdesc)
            rowlist.append(currow)
            i += 1

    neg_n = (maxprobenum - len(rowlist)) / 6
    for neg_idx in range(0,len(validneg)):
        cur_pos_seq = validneg[neg_idx]["sequence"]
        cur_type = validneg[neg_idx]["type"]
        for r in ["r1","r2","r3"]:
            for ori in ["o1","o2"]:
                cur_id = "%s%s" % (arrayname, str(i).zfill(zf))
                curseq = cur_pos_seq if ori == "o1" else bio.revcompstr(cur_pos_seq)
                wprimer = curseq + primer
                rowdesc = "%s_%s_%s" % (cur_type,ori,r)
                currow = "%s\t%s\t%s\tNa|Na\tNa\t%s\tchr1:0-0" % (cur_id, wprimer, rowdesc , rowdesc)
                rowlist.append(currow)
                i += 1
    with open("Coop2Ets.tsv",'w') as f:
        f.write("\n".join(rowlist))

    df_arr = pd.read_csv("Coop2Ets.tsv",header=None,sep="\t")[[1,2]]
    df_arr.columns = ["sequence","id"]
    sanity_check(df_arr, imads12, primer)
    print("done")
