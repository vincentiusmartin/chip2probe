import pandas as pd
import numpy as np

def revcompstr(seq):
    rev = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([rev[base] for base in reversed(seq)])

def count_sim_strs(a,b):
    u=zip(a,b)
    x = 0
    for i,j in u:
        if i==j:
            x += 1
    return x

def make_coopfile(filepath, id_prefix = "Coop1Ets", seqlen = 36, numrep = 3, log=False):
    df = pd.read_csv(filepath, sep="\t")
    dpref = df[df["ID"].str.startswith(id_prefix).fillna(False)][['ID','Name','Sequence','Alexa488Adjusted']]

    if log:
        dpref['Alexa488Adjusted'] = np.log(dpref['Alexa488Adjusted'])

    dpref["Group"] = dpref.apply (lambda row: "_".join(row["Name"].split("_")[:-3]), axis=1)
    g = dpref.groupby("Group")

    reslist = []
    nctrl_list = []
    keyid = {}
    i = 0
    for name, group  in g:
        if i % 100 == 0:
            print("Progress: %d/%d %.2f%%" % (i,len(g),float(i)/len(g) * 100))
        i += 1
        if i == 10:
            break
        for idx, row in group.iterrows():
            cur_name, mut_type, rep, ori = row["Name"].rsplit('_', 3)
            # Checking here for duplicated entry that's why we group per sequence
            curseq = row["Sequence"][:seqlen]
            revcurseq = revcompstr(curseq)

            curkey = ""
            for key in gdict.keys():
                if count_sim_strs(curseq,key) > 30 or count_sim_strs(revcurseq,key) > 30:
                    curkey = key
                    break
            if not curkey:
                curkey = curseq
                keyid[curkey] = len(gdict)
                gdict[curkey] = {}
            if mut_type not in gdict[curkey]:
                if keyid[curkey] > 0:
                    gdict[curkey][mut_type] = {"Name": "%s_dup%d_%s" % (cur_name,keyid[curkey], mut_type)}
                else:
                    gdict[curkey][mut_type] = {"Name": "%s_%s" % (cur_name, mut_type)}
            if "Sequence" not in gdict[curkey] and ori == "o1":
                # we save the positive strand here
                gdict[curkey][mut_type]["Sequence"] = row["Sequence"][:seqlen]
            gdict[curkey][mut_type]["%s_%s" % (ori,rep)] = row["Alexa488Adjusted"]
        for s in gdict:
            for gelm in gdict[s]:
                o1_intensity = [gdict[s][gelm]["o1_r%d" % i] for i in range(1,numrep+1)]
                o2_intensity = [gdict[s][gelm]["o2_r%d" % i]  for i in range(1,numrep+1)]
                gdict[s][gelm]["Median"] = np.median(o1_intensity + o2_intensity)
                gdict[s][gelm]["Median_o1"] = np.median(o1_intensity)
                gdict[s][gelm]["Median_o2"] = np.median(o2_intensity)

        for k in gdict.keys():
            for mtype in gdict[k]:
                if "NegativeCtrl" in name:
                    nctrl_list.append(gdict[k][mtype])
                else:
                    reslist.append(gdict[k][mtype])

    o_r_list =  ["%s_r%s" % (o,r) for o in ["o1","o2"] for r in range(1,numrep + 1)]
    res_df = pd.DataFrame(reslist, columns=["Name","Median","Median_o1","Median_o2"] + o_r_list + ["Sequence"]).sort_values(by=['Name'])
    nctrl_df = pd.DataFrame(nctrl_list, columns=["Name","Median","Median_o1","Median_o2"] + o_r_list + ["Sequence"]).sort_values(by=['Name'])
    return res_df,nctrl_df


if __name__ == "__main__":
    r,n = make_coopfile("ets1_alldata.txt")
    r.to_csv("test.tsv",index=False,sep="\t")
    n.to_csv("negtest.tsv",index=False,sep="\t")
