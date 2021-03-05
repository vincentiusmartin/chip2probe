import pandas as pd
import chip2probe.training_gen.arranalysis as arr
import chip2probe.util.bio as bio
import chip2probe.util.stats_r as st

def check_val(indf, wt, m1, m2, m3):
    df = indf[(indf["Sequence"] == wt) | (indf["Sequence"] == m1) | (indf["Sequence"] == m2) \
            | (indf["Sequence"] == m3) | (indf["Sequence"] == bio.revcompstr(wt)) \
            | (indf["Sequence"] == bio.revcompstr(m1)) | (indf["Sequence"] == bio.revcompstr(m2)) \
            | (indf["Sequence"] == bio.revcompstr(m3))
             ] \
         .sort_values(["ori", "type", "rep"])
    twodict = df[df["type"] == "wt"].groupby('ori')['Alexa488Adjusted'].apply(list).to_dict()
    print(twodict)

    m1df = df[df["type"] == "m1"][['ori','Alexa488Adjusted']]
    m2df = df[df["type"] == "m2"][['ori','Alexa488Adjusted']]
    m3df = df[df["type"] == "m3"][['ori','Alexa488Adjusted']]

    onedf = m1df.merge(m2df, on='ori').merge(m3df, on='ori')
    onedf['indiv'] = onedf["Alexa488Adjusted_x"] + onedf["Alexa488Adjusted_y"] - onedf["Alexa488Adjusted"]
    onedict = onedf.groupby('ori')['indiv'].apply(list).to_dict()

    for ori in ['o1','o2']:
        p = st.wilcox(twodict[ori], onedict[ori], "greater")
        print(ori,p)

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata"

    orig, negorig = arr.read_chamber_file("%s/191030_coop-PBM_Ets1_v1_2nd/2.processed_gpr/20191004_258614510001_ETS1_550_5_1-4_alldata.txt"%basepath, "Coop1Ets", seqcols=["Name","type","rep","ori"], negcols=["Name","rep","ori"])
    cust, negcust = arr.read_chamber_file("%s/210102_validation_array_ets1_v2_2/30nMEts1_alexa488_550_10_alldata.txt"%basepath, "Coop2Ets")

    wt = "TTGCCGGATGAGGTTGGTTGTCTCTAATTTCCTCTC"
    m1 = "TTGCCGCATGAGGTTGGTTGTCTCTAATTTCCTCTC"
    m2 = "TTGCCGGATGAGGTTGGTTGTCTCTAATTGCCTCTC"
    m3 = "TTGCCGCATGAGGTTGGTTGTCTCTAATTGCCTCTC"

    check_val(orig, wt, m1, m2, m3)
