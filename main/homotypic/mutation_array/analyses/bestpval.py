import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import pickle
import numpy as np
import itertools

def get_wtdf(path, sl):
    selected = pd.read_csv(path)
    wtseqs = selected[selected["comment"] == "wt"][["sequence"]].drop_duplicates().rename({"sequence":"Sequence"},axis=1)
    wtseqs = wtseqs.merge(sl).rename({"label":"wtlabel"},axis=1).drop_duplicates().reset_index()
    wtseqs["id"] = wtseqs.index + 1
    wtseqs = wtseqs.sort_values("id")
    wtseqs["id_numeric"] = wtseqs["id"]
    wtseqs["id"] = wtseqs.apply(lambda x: "%s_%s" % (x["id"], x["wtlabel"]), axis=1)
    wtseqs = wtseqs[wtseqs['wtlabel'] != "fail_cutoff"]
    return wtseqs[["Sequence", "wtlabel"]]

def label_row(p1,p2,thres1,thres2):
    if (p1 < thres1 and p2 < thres2) or (p1 < thres2 and p2 < thres1):
        return "cooperative"
    else:
        return "additive"

def get_bestpvals(lbled, wtseqs):
    fl = np.linspace(0.01, 0.5, 20)
    lbled_both = lbled["o1"][["Name", "p"]].merge(lbled["o2"][["Name","p"]], on="Name", suffixes=("_o1", "_o2"))
    check = lbled_both.merge(wtseqs, on="Name")
    cmbs = list(itertools.combinations(fl, 2))
    for c in cmbs:
        check["label"] = check.apply(lambda x: label_row(x["p_o1"],x["p_o2"],c[0],c[1]), axis=1)
        addcount = check[check["wtlabel"] == "additive"].shape[0]
        addtot = check[(check["wtlabel"] == "additive") & (check["label"] == "additive")].shape[0]
        coopcount = check[check["wtlabel"] == "cooperative"].shape[0]
        cooptot = check[(check["wtlabel"] == "cooperative") & (check["label"] == "cooperative")].shape[0]
        alltot = check[(check["wtlabel"] == check["label"])].shape[0]
        percentot = float(alltot)/check.shape[0] * 100
        print("(%.2f,%.2f): additive %d/%d, cooperative %d/%d, total %d/%d (%.2f%%)" % (c[0],c[1],addtot,addcount,cooptot,coopcount,alltot,check.shape[0],percentot))

def get_matching_wts(lbled, wtseqs, pvals):
    print(wtseqs)
    lbled_both = lbled["o1"][["Name", "p"]].merge(lbled["o2"][["Name","p"]], on="Name", suffixes=("_o1", "_o2"))
    check = lbled_both.merge(wtseqs, on="Name")
    check["label"] = check.apply(lambda x: label_row(x["p_o1"],x["p_o2"],pvals[0],pvals[1]), axis=1)
    check_match = check[(check["wtlabel"] == check["label"])]
    return check_match[["Name","Sequence","label"]]


if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/"

    seqlbled = pd.read_csv("%s/chip2probe/output/homotypic/training/seqlbled.csv" % basepath)
    wtdf = get_wtdf("%s/chip2probe/output/array_design_files/Coop2Ets_validation/custom_probes_selected.csv" % basepath, seqlbled)

    # 201128_validation_array_ets1_v2_1/10nMEts1_alexa488_550_20_alldata.txt
    # 210102_validation_array_ets1_v2_2/20nMEts1_alexa488_550_10_alldata.txt
    df, neg = arr.read_chamber_file("%s/probedata/210102_validation_array_ets1_v2_2/30nMEts1_alexa488_550_10_alldata.txt" % basepath, "Coop2Ets")
    wtnames = df.merge(wtdf,on="Sequence")[["Name"]].drop_duplicates()
    wtseqs = df.merge(wtdf,on="Sequence")[["Name","Sequence","wtlabel"]].drop_duplicates().sort_values("Name")
    print(wtseqs["wtlabel"].value_counts())
    cutoff = float(neg.groupby(["Name","ori"]).median().reset_index()[["Alexa488Adjusted"]].quantile(0.95))# confirm
    indiv = pickle.load( open( "%s/chip2probe/output/homotypic/custom_sequences/20nM/indivsum.p" % basepath, "rb" ) )
    two = pickle.load( open( "%s/chip2probe/output/homotypic/custom_sequences/20nM/twosites.p" % basepath, "rb" ) )

    # lbled = arr.label_replicas_permutation(indiv, two, df, cutoff=cutoff, fdrcor=False, pcut=0.05, affcol="Alexa488Adjusted")
    # pickle.dump( lbled, open( "labeled.p", "wb" ) )
    lbled = pickle.load( open( "labeled.p", "rb" ) )
    lbled["o1"] = lbled["o1"].merge(wtnames)
    lbled["o2"] = lbled["o2"].merge(wtnames)

    #get_bestpvals(lbled, wtseqs)

    bestp = [0.11, 0.47]
    cm = get_matching_wts(lbled, wtseqs, bestp)
    cm.to_csv("matchingwts.csv", index=False)
