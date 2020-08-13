import os
os.chdir("../..")

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pbmescore import PBMEscore
from chip2probe.sitespredict.dnasequence import DNASequence

if __name__ == "__main__":
    imads12_paths = ["input/site_models/imads_models/Ets1_w12_GGAA.model", "input/site_models/imads_models/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19) # 0.2128

    escore = PBMEscore("input/site_models/escores/Ets1_8mers_11111111.txt")

    seq = DNASequence("CAGCTGGCCGGAACCTGCGTCCCCTTCCCCCGCCGC", imads12, escore)
    print(seq.sites)
