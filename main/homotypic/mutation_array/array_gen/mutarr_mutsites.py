import os
os.chdir("../../../..")
import pandas as pd
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pbmescore import PBMEscore
from chip2probe.sitespredict.dnasequence import DNASequence
from chip2probe.sitespredict.sitesplotter import SitesPlotter

if __name__ == "__main__":
    df = pd.read_csv("custom_probes_selected.csv")

    # Get array details
    custom_count = df.groupby(["select","muttype"])["seqid"] \
        .count()
    print("----------Number of sequences per group----------")
    print(custom_count)

    wt_count = df.loc[df["comment"] == "wt"]\
        .groupby(["select","muttype"])["seqid"] \
        .count()
    print("----------Number of wt per group----------")
    print(wt_count)

    # Load escore object
    escore = PBMEscore("input/site_models/escores/Ets1_8mers_11111111.txt")

    # Load imads object
    imads12_paths = ["input/site_models/imads_models/Ets1_w12_GGAA.model", "input/site_models/imads_models/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19) # 0.2128

    group_sample = df.drop_duplicates(subset=['sequence']).groupby(["select","muttype"]) \
        .apply(lambda x: x.sample(10)).reset_index(drop=True) \
        .sort_values(by=["select","muttype","seqid"])
    group_sample["type"] = group_sample.apply(lambda x: "%s_%s_%s" % (x["seqid"], x["select"], x["muttype"]),axis=1)
    gs = group_sample[["type","sequence"]].to_dict(orient="records")

    # Get unique samples from each group
    mutplots = {}
    for i in range(len(gs)):
        ds = DNASequence(gs[i]["sequence"], imads12, escore)
        res = {}
        res["wt"] = {"sequence":gs[i]["sequence"],"plt":[]}
        res["m1"] = ds.abolish_sites([0],"to_eliminate",barrier=2,seqonly=False)
        res["m2"] = ds.abolish_sites([1],"to_eliminate",barrier=2,seqonly=False)
        res["m3"] = ds.abolish_sites([0,1],"to_eliminate",barrier=2,seqonly=False)
        newseqs = [res[r]["sequence"] for r in ["wt","m1","m2","m3"]]
        mutpred = imads12.predict_sequences(newseqs,only_pred = True)
        mutlen = [len(mutpred[k]) for k in ['1','2','3','4']]
        if mutlen == [2,1,1,0]: # check that we are not creating new site
            mutplots["%s_wt"%gs[i]["type"]] = res["wt"]
            mutplots["%s_m1"%gs[i]["type"]] = res["m1"]
            mutplots["%s_m2"%gs[i]["type"]] = res["m2"]
            mutplots["%s_m3"%gs[i]["type"]] = res["m3"]

    allseqs = {id: mutplots[id]["sequence"] for id in mutplots}

    escore_pred_list = escore.predict_sequences(allseqs)
    imads_pred_list = imads12.predict_sequences(allseqs, use_threshold=False)

    # Make the plot objects, make_plot_data accepts prediction result
    escore_plot = escore.make_plot_data(escore_pred_list)
    imads_plot = imads12.make_plot_data(imads_pred_list)

    # Generate sequence plot
    sp = SitesPlotter()
    #sp.plot_seq_combine([imads_plot,escore_plot, mutplots], filepath="plot.pdf")

    # --------- make wt m1 m2 m3 for all sequences ---------
    allseqlist = list(df["sequence"].unique())
    print("Sequence count",len(allseqlist))
    allarrseqs = []
    for i in range(len(allseqlist)):
        if i % 100 == 0:
            print("Progress %d/%d" % (i+1,len(allseqlist)))
        ds = DNASequence(allseqlist[i], imads12, escore)
        res = {}
        res["wt"] = allseqlist[i]
        res["m1"] = ds.abolish_sites([0],"to_eliminate",barrier=2,seqonly=True)
        res["m2"] = ds.abolish_sites([1],"to_eliminate",barrier=2,seqonly=True)
        res["m3"] = ds.abolish_sites([0,1],"to_eliminate",barrier=2,seqonly=True)
        newseqs = [res[r] for r in ["wt","m1","m2","m3"]]
        mutpred = imads12.predict_sequences(newseqs,only_pred = True)
        mutlen = [len(mutpred[k]) for k in ['1','2','3','4']]
        if mutlen == [2,1,1,0]: # check that we are not creating new site
            allarrseqs.append(["%d_wt"%(i+1),res["wt"]])
            allarrseqs.append(["%d_m1"%(i+1),res["m1"]])
            allarrseqs.append(["%d_m2"%(i+1),res["m2"]])
            allarrseqs.append(["%d_m3"%(i+1),res["m3"]])

    print("Number of row %d" % len(allarrseqs))

    arrayseqs_df = pd.DataFrame(allarrseqs,columns=["id","sequence"])
    wtcount = arrayseqs_df.loc[arrayseqs_df["id"].str.contains("wt")].shape[0]
    print("Total sequences %d/%d" % (wtcount,len(allseqlist)))
    arrayseqs_df.to_csv("arrayseqs.csv",index=False, header=True)
