'''
Created on Aug 21, 2019

@author: vincentiusmartin
'''
from sitesfinder.pbmescore import PBMEscore
import pickle
from sitesfinder.imadsmodel import iMADSModel
from sitesfinder.imads import iMADS

if __name__ == '__main__':
    escore_short_path = "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/ets1_escores.txt"
    escore_map_path = "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/index_short_to_long.csv"
    escore = PBMEscore(escore_short_path, escore_map_path)
    
    modelcores = ["GGAA",  "GGAT"]
    modelpaths = ["/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_preds/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
                 "/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_preds/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"]
    models = [iMADSModel(modelpath, modelcore, 20, [1,2,3]) for modelpath, modelcore in zip(modelpaths, modelcores)]
    ims = iMADS(models, 0.2128) # 0.2128 is for the ETS1 cutoff

    
    seq = "AAGGCGACGGTTTCCGGTTAGTGGAATCACGGTCCC"
    
    p = escore.predict_sequence(seq)
    i = ims.predict_sequence(seq)
    
    print(i)