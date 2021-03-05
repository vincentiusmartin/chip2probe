import os
os.chdir("../../../..")

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.sitesplotter import SitesPlotter

import chip2probe.training_gen.arranalysis as arr

if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/210102_validation_array_ets1_v2_2"
    df, neg = arr.read_chamber_file("%s/30nMEts1_alexa488_550_10_alldata.txt" % basepath, "Coop2Ets") # 30nMEts1_alexa488_550_10_alldata.txt
    df = df.sort_values(["Name", "ori", "type", "rep"])[["Name", "Sequence", "type", "ori"]].drop_duplicates()
    df["Name"] = df.apply(lambda x: "%s_%s_%s" % (x["Name"],x["ori"],x["type"]), axis=1).head(1000)

    imads_paths = ["input/site_models/imads_models/Ets1_w12_GGAA.model", "input/site_models/imads_models/Ets1_w12_GGAT.model"]
    imads_cores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_paths, imads_cores)]
    imads = iMADS(imads_models, 0.19) # 0.2128

    imadslist = imads.predict_sequences(df, key_colname="Name", sequence_colname="Sequence")
    imadsplot = imads.make_plot_data(imadslist)
    sp = SitesPlotter()
    sp.plot_seq_combine([imadsplot], filepath="plot.pdf")
