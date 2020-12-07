import pandas as pd

import chip2probe.training_gen.arranalysis as arr

def read_chfile(path):
    df = pd.read_csv(path, sep="\t")[["Name", "ID","Sequence","Alexa488Adjusted"]]
    df = df[df["ID"].str.contains("Coop2Ets", na=False)]
    df["Sequence"] = df["Sequence"].str[:36]
    negdf = df[df["Name"].str.contains("NegativeCtrl", na=False)]
    negdf[["Name","ori","rep"]] =  negdf["Name"].str.rsplit("_", n = 2, expand = True)
    df = df[~df["Name"].str.contains("NegativeCtrl", na=False)]
    df[["Name","type","ori","rep"]] = df["Name"].str.rsplit("_", n = 3, expand = True)
    return df, negdf

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/201128_validation_array_ets1_v2_1"
    ch1, neg1 = read_chfile("%s/1nMEts1_alexa488_550_20_alldata.txt" % basepath)
    ch2, neg2 = read_chfile("%s/10nMEts1_alexa488_550_20_alldata.txt" % basepath)
    log = False

    arr.plot_chamber_corr(ch1, ch2, path="custom_corr_550_20.png", median=True, extrajoincols=["type","ori"],
                        namecol="Name", log=log, xlab="Ets1 1nM", ylab="Ets1 10nM")
    arr.plot_chamber_corr(neg1, neg2, path="neg_corr_550_20.png", namecol="Name", extrajoincols=["ori"], median=True,
                        log=log, xlab="Ets1 1nM", ylab="Ets1 10nM")

    # get the negative control cutoff, we do it from first chamber
    cutoff = neg2.groupby(["Name","ori"]).median().reset_index()[["Alexa488Adjusted"]]. quantile(0.95)
    print("Negative control cutoff", cutoff)
