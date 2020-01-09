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
import itertools
from cooperative import coopfilter



def main():
    curdir = "/Users/vincentiusmartin/Research/chip2gcPBM/"
    analysis_path = curdir + "result/ets1_HepG2/analysis_result/"
    
    
    with open(analysis_path + "sitefiles_list.txt", 'r') as f:
        sitelist = [line.strip() for line in f.readlines()]
    
    #slist = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/../result/ets1_A549/analysis_result/sites_within_d3_span50.tsv"
    slist = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/../result/cistrome/cistrome_ets1_37927/analysis_result/sites_within_d2_span100.tsv"
    seqdf = pd.read_csv(slist, sep='\t')

    escore_short_path = "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/ets1_escores.txt"
    escore_map_path = "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/index_short_to_long.csv"
    escore = PBMEscore(escore_short_path, escore_map_path)
    es_preds = escore.predict_sequences(seqdf)
    esplots = escore.plot(es_preds)
    


    modelcores = ["GGAA",  "GGAT"]
    modelpaths = ["/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_preds/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
                 "/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_preds/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"]
    models = [iMADSModel(modelpath, modelcore, 20, [1,2,3]) for modelpath, modelcore in zip(modelpaths, modelcores)]
    ims = iMADS(models, 0.2128) # 0.2128 is for the ETS1 cutoff
    ims_preds = ims.predict_sequences(seqdf)
    #pickle.dump(ims_preds, open("ims_preds.pickle","wb"))
    #ims_preds = pickle.load(open("ims_preds.pickle","rb"))
    
    imadsplots = ims.plot(ims_preds)
    
    pc = PlotCombiner()
    pc.plot_seq_combine([imadsplots,esplots], filepath="plot.pdf")
    
    filtered_sites = {}
    print("Site filtering...")
    for key in es_preds:
    #for key in ["sequence12"]:
    # loose filtering in here
        bs = Sequence(es_preds[key],ims_preds[key],escore_cutoff=0.35, escore_gap = 2)
        if bs.site_count() == 2:
            filtered_sites[key] = bs
    #site_list = [{**{"key":site, "sequence":filtered_sites[site].sequence},**filtered_sites[site].get_sites_dict()} for site in filtered_sites]
    #columns = ["key", "site_start_1", "site_start_2", "site_end_1", "site_end_2", "site_pos_1", "site_pos_2", "imads_score_1", "imads_score_2", "sequence"
    
    #pickle.dump(filtered_sites, open("test_fsites2.pickle","wb"))
    #filtered_sites = pickle.load(open("test_fsites2.pickle","rb"))
    ###
    seqdict = {}
    funcdict = {}
    filtered_probes = []
    # TODO: tmr look at 110,271
    for key in filtered_sites:
    #for key in ["sequence12"]:
        # Visualization part
        seqdict["%s-wt" % key] = filtered_sites[key].sequence
        for idx,mut in enumerate([[0],[1],[0,1]]): 
            mutseq = filtered_sites[key].abolish_sites(mut,escore)
            seqdict["%s-m%d" % (key,idx + 1)] = mutseq.sequence
            funcdict["%s-m%d" % (key,idx + 1)] = mutseq.plot_functions
        for e in list(itertools.product([0,1,2],[0.41, 0.4, 0.35])): # ,0.35
            egapthres = e[0]
            ecutoff = e[1]
            if coopfilter.filter_coopseq(seqdict["%s-wt"%key], seqdict["%s-m1"%key],
                                     seqdict["%s-m2"%key], seqdict["%s-m3"%key],
                                     filtered_sites[key].get_sites_dict(), escore,
                                     escore_cutoff=ecutoff, escore_gap = egapthres):
                filtered_probes.append({"key":key, "wt":seqdict["%s-wt"%key], "m1":seqdict["%s-m1"%key],
                                    "m2":seqdict["%s-m2"%key], "m3":seqdict["%s-m3"%key], "ecutoff": ecutoff, 
                                    "egapthres":egapthres,"distance":filtered_sites[key].get_sites_dist()})
                break

    pp = escore.plot(escore.predict_sequences(seqdict),additional_functions=funcdict)
    pc.plot_seq_combine([pp], filepath="plot-mut.pdf")
    
    if(filtered_probes):
        pd.DataFrame(filtered_probes).to_csv("mutated_probes.tsv",sep="\t",index=False,columns=["key","wt","m1","m2","m3","ecutoff", "egapthres","distance"])
    
    print("Done")

if __name__ == '__main__':
    main()
    