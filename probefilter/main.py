'''
Created on Jul 16, 2019

@author: vincentiusmartin
'''

import pandas as pd

from sitesfinder.imads import iMADS
from sitesfinder.imadsmodel import iMADSModel
from sitesfinder.plotcombiner import PlotCombiner
from sitesfinder.pbmescore import PBMEscore
from sitesfinder.sequence import Sequence
from sitesfinder.prediction.basepred import BasePrediction

import pickle


def main():
    curdir = "/Users/vincentiusmartin/Research/chip2gcPBM/"
    analysis_path = curdir + "result/ets1_HepG2/analysis_result/"
    
    
    with open(analysis_path + "sitefiles_list.txt", 'r') as f:
        sitelist = [line.strip() for line in f.readlines()]
    
    seqdf = pd.read_csv(sitelist[0], sep='\t')
    
    escore_short_path = "/Users/vincentiusmartin/Research/chip2gcPBM/escores/ets1_escores.txt"
    escore_map_path = "/Users/vincentiusmartin/Research/chip2gcPBM/escores/index_short_to_long.csv"
    escore = PBMEscore(escore_short_path, escore_map_path)
    es_preds = escore.predict_sequences(seqdf)
    #esplots = escore.plot(es_preds)
    """
    modelcores = ["GGAA",  "GGAT"]
    modelpaths = ["/Users/vincentiusmartin/Research/chip2gcPBM/imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
                 "/Users/vincentiusmartin/Research/chip2gcPBM/imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"]
    models = [iMADSModel(modelpath, modelcore, 20, [1,2,3]) for modelpath, modelcore in zip(modelpaths, modelcores)]
    ims = iMADS(models, 0.2128) # 0.2128 is for the ETS1 cutoff
    ims_preds = ims.predict_sequences(seqdf)
    #imadsplots = ims.plot(ims_preds)
    
    pc = PlotCombiner()
    #pc.plot_seq_combine([imadsplots,esplots], filepath="plot.pdf")
    
    filtered_sites = {}
    print("Site filtering...")
    for key in es_preds:
        es_pred1 = es_preds[key]
        ims_pred1 =  ims_preds[key]
        bs = Sequence(ims_pred1,es_pred1,es_preds[key].sequence)
        if bs.site_count() == 2:
            filtered_sites[key] = bs
    site_list = [{**{"key":site, "sequence":filtered_sites[site].sequence},**filtered_sites[site].get_sites_dict()} for site in filtered_sites]
    columns = ["key", "site_start_1", "site_start_2", "site_end_1", "site_end_2", "site_pos_1", "site_pos_2", "imads_score_1", "imads_score_2", "sequence"]
    pd.DataFrame(site_list).to_csv("out.csv", index=False, columns=columns, float_format='%.4f')
    """
    
    #pickle.dump(filtered_sites, open("test_fsites.pickle","wb"))
    filtered_sites = pickle.load(open("test_fsites.pickle","rb"))
    pc = PlotCombiner()
    ###
    seqdict = {}
    funcdict = {}
    for key in filtered_sites:
        seqdict["%s-wild" % key] = filtered_sites[key].sequence
        for idx,mut in enumerate([[0],[1],[0,1]]): #  # [0],[1],
            mutseq = filtered_sites[key].mutate_sites(mut,escore)
            seqdict["%s-mut%d" % (key,idx + 1)] = mutseq.sequence
            funcdict["%s-mut%d" % (key,idx + 1)] = mutseq.plot_functions

    pp = escore.plot(escore.predict_sequences(seqdict),additional_functions=funcdict)
    pc.plot_seq_combine([pp], filepath="plot-mut.pdf")
    
    print("Done")
    

if __name__ == '__main__':
    main()
    