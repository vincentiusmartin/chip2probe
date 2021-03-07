import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import numpy as np
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.sitesplotter import SitesPlotter

def get_wtdf(path, sl):
    selected = pd.read_csv(path)
    wtseqs = selected[selected["comment"] == "wt"][["sequence"]].drop_duplicates() #.rename({"sequence":"Sequence"},axis=1)
    wtseqs = wtseqs.merge(sl).rename({"label":"wtlabel"},axis=1).drop_duplicates().reset_index()
    wtseqs["id"] = wtseqs.index + 1
    wtseqs = wtseqs.sort_values("id")
    wtseqs["id_numeric"] = wtseqs["id"]
    wtseqs["id"] = wtseqs.apply(lambda x: "%s_%s" % (x["id"], x["wtlabel"]), axis=1)
    wtseqs = wtseqs[wtseqs['wtlabel'] != "fail_cutoff"]
    return wtseqs

def make_wtseq_col(df):
    onlywt = df[(df["type"] == "wt") & (df["ori"] == "o1")][["Name","Sequence"]].drop_duplicates()
    df = df.drop(["Sequence"], axis=1)
    df = df.merge(onlywt, on="Name")
    df = df.sort_values(["Name","type","ori","rep"])
    df["Alexa488Adjusted"] = np.log(df["Alexa488Adjusted"])
    return df

def intersectwt(path, wtdf):
    df = pd.read_csv(path)
    return wtdf.merge(df, on="Sequence")

def get_wtmt(indf, wtdf):
    seqid = indf.merge(wtdf, on="Sequence")[["Name", "id_numeric"]].sort_values("id_numeric").drop_duplicates()
    dfwtmt = seqid.merge(indf, on="Name")
    dfwtmt = dfwtmt[dfwtmt["ori"] == "o1"][["Sequence","type","id_numeric"]].drop_duplicates()
    dfwtmt["type"] = pd.Categorical(dfwtmt['type'], ["wt", "m1", "m2", "m3"])
    dfwtmt = dfwtmt.sort_values(["id_numeric","type"])
    dfwtmt = dfwtmt.groupby("id_numeric").filter(lambda x: len(x) <= 4)
    dfwtmt["id_numeric"] = dfwtmt.apply(lambda x: "%s_%s" % (x["id_numeric"], x["type"]), axis=1)
    return dfwtmt

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/"

    seqlbled = pd.read_csv("%s/chip2probe/output/homotypic/training/seqlbled.csv" % basepath)
    wtdf = get_wtdf("%s/chip2probe/output/array_design_files/Coop2Ets_validation/custom_probes_selected.csv" % basepath, seqlbled)

    origdf, neg_orig = arr.read_chamber_file("%s/probedata/191030_coop-PBM_Ets1_v1_2nd/2.processed_gpr/20191004_258614510001_ETS1_550_5_1-4_alldata.txt"%basepath, seqcols=["Name","type","rep","ori"], negcols=["Name","rep","ori"], key="Coop1Ets")
    origdf[["Sequence","type","ori"]].drop_duplicates().to_csv("seqsorig.csv",index=False)
    import sys
    sys.exit()
    cust10df, neg10_cust = arr.read_chamber_file("%s/probedata/201128_validation_array_ets1_v2_1/10nMEts1_alexa488_550_20_alldata.txt"%basepath, key="Coop2Ets")
    cust20df, neg20_cust = arr.read_chamber_file("%s/probedata/210102_validation_array_ets1_v2_2/20nMEts1_alexa488_550_10_alldata.txt"%basepath, key="Coop2Ets")
    cust30df, neg30_cust = arr.read_chamber_file("%s/probedata/210102_validation_array_ets1_v2_2/30nMEts1_alexa488_550_10_alldata.txt"%basepath, key="Coop2Ets")

    imads_paths = ["%s/chip2probe/input/site_models/imads_models/Ets1_w12_GGAA.model" % basepath,
                    "%s/chip2probe/input/site_models/imads_models/Ets1_w12_GGAT.model" % basepath]
    imads_cores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_paths, imads_cores)]
    imads = iMADS(imads_models, 0.19) # 0.2128

    orig_pred = imads.predict_sequences(get_wtmt(origdf, wtdf), key_colname="id_numeric", sequence_colname="Sequence")
    orig_plot = imads.make_plot_data(orig_pred)

    cust_pred = imads.predict_sequences(get_wtmt(cust10df, wtdf), key_colname="id_numeric", sequence_colname="Sequence")
    cust_plot = imads.make_plot_data(cust_pred)

    sp = SitesPlotter()
    sp.plot_seq_combine([orig_plot], filepath="origplot.pdf")
    sp.plot_seq_combine([cust_plot], filepath="custplot.pdf")

    import sys
    sys.exit()

    origdf = wtdf.merge(make_wtseq_col(origdf), on="Sequence")
    origdf["arrtype"] = "orig"
    cust10df = wtdf.merge(make_wtseq_col(cust10df), on="Sequence")
    cust10df["arrtype"] = "10nM"
    cust10df[["Sequence","id"]].to_csv("custseqswt.csv",index=False)
    cust20df = wtdf.merge(make_wtseq_col(cust20df), on="Sequence")
    cust20df["arrtype"] = "20nM"
    cust30df = wtdf.merge(make_wtseq_col(cust30df), on="Sequence")
    cust30df["arrtype"] = "30nM"

    alldf = pd.concat([origdf, cust10df, cust20df, cust30df])
    alldf["arrtype"] = pd.Categorical(alldf['arrtype'], ["orig", "10nM", "20nM", "30nM"])
    alldf = alldf.sort_values(["id_numeric","arrtype"])
    alldf["id"] = alldf.apply(lambda x: "%s_%s" % (x["id"], x["arrtype"]), axis=1)

    arr.plot_multi_scatterbox("cust_alltypes.pdf", alldf, ["wt","m1","m2","m3"], namecol="id", affcol="Alexa488Adjusted", pline=False, allplot=False)

    ind_orig = intersectwt("%s/chip2probe/output/homotypic/training/indiv.csv" % basepath, wtdf)
    ind_orig["arrtype"] = "orig"
    two_orig = intersectwt("%s/chip2probe/output/homotypic/training/two.csv" % basepath, wtdf)
    two_orig["arrtype"] = "orig"

    ind10_cust = intersectwt("%s/chip2probe/output/homotypic/custom_sequences/10nM/indiv.csv" % basepath, wtdf)
    ind10_cust["arrtype"] = "10nM"
    two10_cust = intersectwt("%s/chip2probe/output/homotypic/custom_sequences/10nM/two.csv" % basepath, wtdf)
    two10_cust["arrtype"] = "10nM"

    ind20_cust = intersectwt("%s/chip2probe/output/homotypic/custom_sequences/20nM/indiv.csv" % basepath, wtdf)
    ind20_cust["arrtype"] = "20nM"
    two20_cust = intersectwt("%s/chip2probe/output/homotypic/custom_sequences/20nM/two.csv" % basepath, wtdf)
    two20_cust["arrtype"] = "20nM"

    ind30_cust = intersectwt("%s/chip2probe/output/homotypic/custom_sequences/30nM/indiv.csv" % basepath, wtdf)
    ind30_cust["arrtype"] = "30nM"
    two30_cust = intersectwt("%s/chip2probe/output/homotypic/custom_sequences/30nM/two.csv" % basepath, wtdf)
    two30_cust["arrtype"] = "30nM"

    allind = pd.concat([ind_orig, ind10_cust, ind20_cust, ind30_cust])
    allind["arrtype"] = pd.Categorical(allind['arrtype'], ["orig", "10nM", "20nM", "30nM"])
    allind = allind.sort_values(["id_numeric","arrtype"])
    allind["id"] = allind.apply(lambda x: "%s_%s" % (x["id"], x["arrtype"]), axis=1)

    alltwo = pd.concat([two_orig, two10_cust, two20_cust, two30_cust])
    alltwo["arrtype"] = pd.Categorical(alltwo['arrtype'], ["orig", "10nM", "20nM", "30nM"])
    alltwo = alltwo.sort_values(["id_numeric","arrtype"])
    alltwo["id"] = alltwo.apply(lambda x: "%s_%s" % (x["id"], x["arrtype"]), axis=1)

    arr.plot_ori_inconsistency(allind, alltwo, namecol="id", prefix_path="all", log=True, allplot=False)
