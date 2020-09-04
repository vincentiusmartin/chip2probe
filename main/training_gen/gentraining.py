import os
os.chdir("../..")
import pandas as pd

from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

def get_sites(seq):
    return seq

if __name__ == "__main__":
    df = pd.read_csv("output/heterotypic/EtsRunx_v1/sequence_labeled_normalized.tsv", sep="\t")

    smdir = "input/site_models/imads_models"

    imads_ets_paths = ["%s/Ets1_w12_GGAA.model" % smdir, "%s/Ets1_w12_GGAT.model" % smdir]
    imads_ets_cores = ["GGAA", "GGAT"]
    imads_ets_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_ets_paths, imads_ets_cores)]
    imads_ets = iMADS(imads_ets_models, 0) # 0.2128

    imads_runx_paths = ["%s/Runx1_w12_GAGGT.model" % smdir, "%s/Runx1_w12_GCGGC.model" % smdir, "%s/Runx1_w12_GCGGT.model" % smdir, "%s/Runx1_w12_GTGGC.model" % smdir, "%s/Runx1_w12_GTGGT.model" % smdir]
    imads_runx_cores = ["GAGGT", "GCGGC", "GCGGT", "GTGGC", "GTGGT"]
    imads_runx_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_runx_paths, imads_runx_cores)]
    imads_runx = iMADS(imads_runx_models, 0.3061) # 0.2128

    kompas_ets = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)
    kompas_runx = Kompas("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/input/site_models/kompas/Runx1_kmer_alignment.txt",
                    core_start = 12, core_end = 17, core_center = 14)

    # df = df["Sequence"]lambda x: get_sites(x))
    # print(df)

    s = "TCTTCCTCCTCGCGGTCGCGGCCGGACGGAGGGTGG" #"CTAGCAACCGCGGGAAGGGGGCGTGGCCGGAAGCCT"
    print("imads runx",imads_runx.predict_sequence(s))
    print("imads ets",imads_ets.predict_sequence(s))

    print("kompas runx",kompas_runx.predict_sequence(s))
    print("kompas ets",kompas_ets.predict_sequence(s))
