import pandas as pd

from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.kompas import Kompas
import chip2probe.training_gen.traingen as tg
from chip2probe.util import bio

def processdf(path):
    df = pd.read_csv(path, sep="\t", header=None)
    newdf = pd.DataFrame()
    newdf["Sequence"] = df[1]
    newdf["type"] = df[2]
    newdf["id"] = df[2].apply(lambda x: x.rsplit("_",3)[0])
    newdf["muttype"] = df[2].apply(lambda x: x.split("_")[2])
    newdf["ori"] = df[2].apply(lambda x: x.split("_")[3])
    return newdf

def join_back(selected, valid_df, map, primerlen, keyid):
    # TODO: check all groups have wt
    vd = valid_df.copy()
    vd = vd[vd["id"].str.contains(keyid)]
    vd["Sequence"] = vd["Sequence"].apply(lambda x: x[:-primerlen])
    vd = vd[(vd["muttype"] == "wt") & (vd["ori"] == "o1")][["Sequence"]].drop_duplicates()
    vd["Sequence"] = vd.replace({"Sequence": map})
    print("Seqlen before joined %d" % vd.shape[0])
    joined = vd.merge(selected, on="Sequence")
    print("Seqlen after joined %d" % joined[["Sequence"]].drop_duplicates().shape[0])
    return joined

def select_probes(df, negdict, numspots = 176024, numrep=3, arrname="Coop3Ets_", primerseq="GTCTTGATTCGCTTGACGCTGCTG"):
    len_d = sum(len(v) for v in negdict.values())
    negspots = len_d * numrep * 2
    availspots = numspots - negspots

    copydf = df.copy()
    copydf[['exp', 'id']] = copydf['id'].str.split('_', 1, expand=True)

    arraylines = []
    id = 0
    maxidlen = 15
    zf = maxidlen - len(arrname)

    exp_group = copydf.groupby("exp")
    exp_group = {e:df_egrp for e,df_egrp in exp_group}
    chunk = availspots // len(exp_group)
    initspot = int(chunk)

    for exp_name in ["er","ee"]:
        exp_df = exp_group[exp_name]
        cur_arrlines = []
        for probe_name, probe_group in exp_df.groupby("id"):
            if len(cur_arrlines) + probe_group.shape[0] > initspot:
                break
            for idx, row in probe_group.iterrows():
                id += 1
                cur_id = "%s%s" % (arrname, str(id).zfill(zf))
                cur_arrlines.append("%s\t%s\t%s\tNa|Na\tNa\t%s\tchr1:0-0" % (cur_id, row["Sequence"], row["type"], row["type"]))
        print("%d spots for %s" % (len(cur_arrlines),exp_name))
        arraylines.extend(cur_arrlines)
        initspot = int(chunk) if availspots - len(arraylines) > 2*initspot else availspots - len(arraylines)

    for k in negdict:
        neglist = negdict[k]
        for i in range(len(neglist)):
            for o in ["o1", "o2"]:
                if o == "o1":
                    seq = str(neglist[i]) + primerseq
                else:
                    seq = bio.revcompstr(neglist[i]) + primerseq
                for r in range(numrep):
                    id += 1
                    cur_id = "%s%s" % (arrname, str(id).zfill(zf))
                    type = "NegCtrl-%s_%d_%s_%s_r%d" % (k, i+1, "neg", o, r)
                    arraylines.append("%s\t%s\t%s\tNa|Na\tNa\t%s\tchr1:0-0" % (cur_id, seq, type, type))

    with open("finarr.txt",'w') as f:
        f.write("\n".join(arraylines))

def validate_array(dfall, kompas_ets, kompas_runx, pwm_ets, pwm_runx, primer="GTCTTGATTCGCTTGACGCTGCTG"):
    dfsub = dfall.drop_duplicates("Sequence").reset_index(drop=True)
    dfsub_groups = dfsub.groupby("id")

    md = len(dfsub_groups) // 100
    i = 0
    valid_ids = []
    for name, group in dfsub_groups:
        if i % md == 0:
            print("Progress validating %d/%d" % (i,len(dfsub_groups)))
        i += 1
        ok = True
        for idx, row in group.iterrows():
            kompdict = {"ets1":kompas_ets}
            pwmdict = {"ets1":pwm_ets}
            if "er_" in row["type"]:
                toelims = [0,1]
                kompdict["runx1"] = kompas_runx
                pwmdict["runx1"] = pwm_runx
            pred = tg.pred_all(row["Sequence"], kompdict, pwmdict)
            numsites = len(pred)
            if  (row["muttype"] == "wt" and numsites != 2) or \
                (row["muttype"] == "m1" and numsites != 1) or \
                (row["muttype"] == "m2" and numsites != 1) or \
                (row["muttype"] == "m3" and numsites != 0) or \
                (row["muttype"] == "neg" and numsites != 0):
                print("Error for sequence %s" % row["Sequence"])
                print(i, numsites, row["muttype"] )
                print(pred)
                ok = False
            elif not row["Sequence"].endswith(primer):
                print("Not ending with primer %s" % row["Sequence"])
                ok = False
            if not ok:
                break
        if ok:
            valid_ids.append({"id":row["id"]})
    valid_ids = pd.DataFrame(valid_ids)
    valid_df = dfall.merge(valid_ids, on='id')
    valid_df.to_csv("valid.csv", index=False)
    return valid_df

def clean_neg(negpath, primer, kompdict, pwmdict, tfname, negcount=-1,):
    df = pd.read_csv(negpath)
    df = df[df["ID"].str.contains("NegCtrl")]
    seqlist = df["Sequence"].tolist()
    print("initneg",len(seqlist))

    # kompdict = {"tf":kompas}
    # pwmdict = {"tf":pwm}

    oklist = []
    nc = len(seqlist) if negcount <= 0 else negcount
    for seq in seqlist:
        if len(oklist) >= nc:
            break
        seqprim = seq + primer
        revprim = bio.revcompstr(seq) + primer
        fwdpred = tg.pred_all(seqprim, kompdict, pwmdict)
        rcpred = tg.pred_all(revprim, kompdict, pwmdict)
        if len(fwdpred) == 0 and len(rcpred) == 0:
            oklist.append(seq)
    print("neglen %s" % tfname,len(oklist))
    return oklist

pd.set_option("display.max_columns",None)
if __name__ == "__main__":

    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"
    pwm_ets = PWM("%s/input/sitemodels/pwm/ets1.txt" % basepath, log=True, reverse=False)
    kompas_ets = Kompas("%s/input/sitemodels/kompas/Ets1_kmer_alignment.txt" % basepath, core_start = 11, core_end = 15, core_center = 12)
    pwm_runx = PWM("%s/input/sitemodels/pwm/runx1.txt" % basepath, 8, 17, log=True, reverse=True)
    kompas_runx = Kompas("%s/input/sitemodels/kompas/Runx1_kmer_alignment.txt" % basepath, core_start = 12, core_end = 17, core_center = 14)

    primer = "GTCTTGATTCGCTTGACGCTGCTG"
    ee_path, er_path = "output/arr_etsets.txt", "output/arr_etsrunx.txt"
    ee_lit_path, er_lit_path = "output/arr_lit_etsets.txt", "output/arr_lit_etsrunx.txt"

    df_ee = processdf(ee_path)
    df_er = processdf(er_path)
    df_lit_ee = processdf(ee_lit_path)
    df_lit_er = processdf(er_lit_path)
    dfall = pd.concat([df_lit_er,df_er,df_lit_ee,df_ee])

    # valid_df = validate_array(dfall, kompas_ets, kompas_runx, pwm_ets, pwm_runx)
    # valid_df.to_csv("output/valid.csv",index=False)
    valid_df = pd.read_csv("output/valid.csv")

    negcount = 253
    negctrl_ets = clean_neg("input/imads_ets1.csv", primer, {"ets1":kompas_ets, "runx1":kompas_runx}, {"ets1":pwm_ets, "runx1":pwm_runx}, "ets", negcount)
    negctrl_runx = clean_neg("input/imads_runx1.csv", primer, {"ets1":kompas_ets, "runx1":kompas_runx}, {"ets1":pwm_ets, "runx1":pwm_runx}, "runx", negcount)
    negctrl = {"ets":negctrl_ets, "runx": negctrl_runx}

    select_probes(valid_df, negctrl, primerseq=primer)

    print("Validating final array...")

    dffin = processdf("finarr.txt")
    # validate_array(dffin, kompas_ets, kompas_runx, pwm_ets, pwm_runx)

    # ======= Join back =======

    selected_ee = pd.read_csv("output/custom_selected_etsets.csv")
    selected_ee_lit =  pd.read_csv("output/lit_sitedata_etsets.csv")
    selected_ee = pd.concat([selected_ee_lit,selected_ee])
    junc_ee = pd.read_csv("output/cleanjunction_map_etsets.csv")
    junc_ee = pd.Series(junc_ee["old"].values,index=junc_ee["new"]).to_dict()

    selected_er = pd.read_csv("output/custom_selected_etsrunx.csv")
    selected_er_lit =  pd.read_csv("output/lit_sitedata_etsrunx.csv")
    selected_er = pd.concat([selected_er_lit,selected_er])
    junc_er = pd.read_csv("output/cleanjunction_map_etsrunx.csv")
    junc_er = pd.Series(junc_er["old"].values,index=junc_er["new"]).to_dict()

    valid_df = pd.read_csv("output/valid.csv")
    df_fin_ee = join_back(selected_ee,dffin,junc_ee,len(primer),"ee").sort_values("id")
    print(df_fin_ee["muttype"].value_counts())
    df_fin_ee.to_csv("selected_ee_fin.csv", index=False)
    df_fin_er = join_back(selected_er,dffin,junc_er,len(primer),"er").sort_values("id")
    print(df_fin_er["muttype"].value_counts())
    df_fin_er.to_csv("selected_er_fin.csv", index=False)

    print("done")
