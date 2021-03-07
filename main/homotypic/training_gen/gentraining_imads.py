from chip2probe.modeler.cooptrain import CoopTrain
import pandas as pd

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"
    # using custom imads model
    imads_paths = ["%s/input/site_models/imads_models/Ets1_w12_GGAA.model" % basepath, "%s/input/site_models/imads_models/Ets1_w12_GGAT.model" % basepath]
    imads_cores = ["GGAA", "GGAT"]
    imads_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_paths, imads_cores)]
    imads = iMADS(imads_models, 0.19) # 0.2128

    df = pd.read_csv("%s/output/homotypic/training/seqlbled.csv" % basepath).drop_duplicates()
    print(df.shape[0])
    df = df[(df["label"] == "cooperative") | (df["label"] == "independent")].rename({"Sequence":"sequence"}, axis=1)
    print(df["label"].value_counts())
    ct = CoopTrain(df["sequence"].tolist(), corelen=4, flip_th=True, imads=imads, ignore_sites_err=True)
    train_df = ct.df.merge(df, on="sequence")
    print(train_df["label"].value_counts())
    print(train_df.shape[0])
    train_df.to_csv("training.csv",index=False)
