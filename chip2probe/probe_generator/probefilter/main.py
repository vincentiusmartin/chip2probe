'''
This is the main file used in running probefilter.

The output is a file of clean non custom probes as well as the optional
clean custom probes, with primers

Authors: Farica Zhuang, Vincentius Martin
Created on Sep 30, 2019
'''

# these paths need to be handled better
import sys
import pandas as pd
import itertools
import os
sys.path.append(os.path.realpath("sitespredict/libsvm/python"))
from chip2probe.probe_generator.probefilter import sitespredict
from sitespredict.pbmescore import PBMEscore
from sitespredict.imadsmodel import iMADSModel
from sitespredict.sitesplotter import SitesPlotter
from sitespredict.imads import iMADS
from sitespredict.sequence import Sequence
from sitespredict.kompas import Kompas
from cooperative import coopfilter
from chip2probe.util import bio as bio

if __name__ == '__main__':
    #-----------------------provided by user-------------------------
    # path to all_sites.tsv
    analysis_path = ""
    # provide either escore_long_path OR (escore_short_path AND escore_map_path)
    # Download escore short and map from Qbic http://qbic.genome.duke.edu/downloads
    # next to 'Download all e-score tables used in our prediction'
    # read escore file
    escore_short_paths = {'ets1': "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/ets1_escores.txt"}
    escore_map_paths = {'ets1': "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/index_short_to_long.csv"}

    # path will be provided by user
    escore_long_paths = {'ets1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/escores/Ets1_8mers_11111111.txt",
                         'runx1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/escores/Runx1_8mers_11111111.txt"}

    #curdir = "/Users/vincentiusmartin/Research/chip2gcPBM/"
    #analysis_path = curdir + "result/cistrome_ets1_37927/analysis_result/"

    # set the model used, either kompas or imads
    method = "imads"
    kmer_align_paths = {'ets1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/kompas/ets1/Ets1_kmer_alignment.txt",
                        'runx1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/kompas/runx1/Runx1_kmer_alignment.txt"}
    modelcores = {'ets1': ["GGAA", "GGAT"],
                  'runx1': ["GAGGT", "GCGGC", "GCGGG", "GCGGT", "GTGGC", "GTGGG", "GTGGT"]}
    modelpaths = {'ets1': ["/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
                           "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"],
                  'runx1': ["/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GAGGT_1a2a3mer_format.model",
                            "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGC_1a2a3mer_format.model",
                            "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGG_1a2a3mer_format.model",
                            "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGT_1a2a3mer_format.model",
                            "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGC_1a2a3mer_format.model",
                            "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGG_1a2a3mer_format.model",
                            "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGT_1a2a3mer_format.model"]
                 }
    imads_cutoff = {'ets1': 0.2128,
                    'runx1': 0.3061}
    proteins = ['ets1']
    # True to indicate mutating sequences to clean them
    mutate = False
    # list of escore gaps
    egaps = [0]
    # list of escore thresholds
    thresholds = [0.4]
    # cutoff for significant escore. If mutate is True, any non-centered sites below this cutoff are
    # mutated and removed
    mutate_cutoff = 0.38
    # gap used in determining what to mutate
    mutate_gap = 0
    # imads would usually use flanks to predict
    predict_flanks = False
    # boolean whether to show flanks used by the model. Kompas has no flanks, only imads does
    show_model_flanks = False
    # user can determine colors for the proteins as a list of tuples (escore color, model color)
    colors = []
    #-------------------------------------------------------------------------------------
    if len(colors) == 0:
        colors = [('crimson', 'plum'), ('steelblue', 'lightblue')]
    # create escore and model objects for each protein
    escores = {}
    models = {}
    for protein in proteins:
        escores[protein] = PBMEscore(escore_long_paths[protein])
        # initialize binding site prediction model used based on method
        if method=='imads':
            # read imads prediction
            imads_models = [iMADSModel(modelpath, modelcore, 20, [1, 2, 3])
                        for modelpath, modelcore in zip(modelpaths[protein], modelcores[protein])]
            models[protein] = iMADS(imads_models, imads_cutoff[protein])  # 0.2128 is for the ETS1 cutoff
        elif method == 'kompas':
            models[protein] = Kompas(protein=protein, threshold=mutate_cutoff,
                                     kmer_align_path=kmer_align_paths[protein])
        # raise exception if model is not supported
        else:
            raise Exception("Supported methods are imads and kompas. Got {}".format(method))

    # read the input dataframe
    pd.set_option('display.max_columns', None)
    sitepath = "%ssites_all.tsv" % analysis_path
    seqdf = pd.read_csv(sitepath, sep='\t')

    # get all combinations of number of sites in peak and peak length as a list of tuples
    spcomb = sorted(list(set(zip(seqdf["sites.in.peak"], seqdf["peaklen"]))), key=lambda x: (x[0],x[1]))

    filtered_probes = []  # we will fill this

    # for each combination, get filtered probes
    for comb in spcomb:
        es_preds = {}
        esplots = {}
        model_preds = {}
        model_plots = {}
        sitenum = comb[0]
        peaklen = comb[1]
        print("Working on sitenum %d, peaklen %d" % (sitenum,peaklen))
        # get sequences with the specified sites in peak and peak length
        curdf = seqdf.loc[(seqdf["sites.in.peak"] == sitenum) & (seqdf["peaklen"] == peaklen)]
        curdf = curdf[:20]
        # get escore and model predictions for each protein
        for protein in proteins:
            protein_num = proteins.index(protein)
            es_preds[protein] = escores[protein].predict_sequences(curdf, key_colname="key")
            esplots[protein] = escores[protein].make_plot_data(es_preds[protein], color=colors[protein_num][0])

            model_preds[protein] = models[protein].predict_sequences(curdf, 
                                                                     key_colname="key", 
                                                                     predict_flanks=predict_flanks)
            model_plots[protein] = models[protein].make_plot_data(model_preds[protein], 
                                                                  color=colors[protein_num][1],
                                                                  show_model_flanks=show_model_flanks)

        # Generate plots
        sp = SitesPlotter()
        # if need to plot, uncomment this
        sp.plot_seq_combine([esplots, model_plots],
                            filepath="%ssitesplot_d%d_p%d.pdf" %
                            (analysis_path, sitenum, peaklen))

        # get filtered sequences
        filtered_seqs = {}
        flanks = {}
        print("Site filtering...")
        print("Number of sites before mutating:", len(es_preds[proteins[0]]))

        # get sequences with 2 significant binding sites
        sites_mutated = 0
        sites_removed = 0
        failed_mutations = 0
        for key in es_preds[proteins[0]]:
            curr_es_preds = {}
            curr_model_preds = {}
            for protein in proteins:
                curr_es_preds[protein] = es_preds[protein][key]
                curr_model_preds[protein] = model_preds[protein][key]
            bs = Sequence(curr_es_preds, curr_model_preds, proteins=proteins,
                          escore_cutoff=mutate_cutoff, escore_gap=mutate_gap,
                          pbmescore=escores)
            if bs.is_valid(mutate):
                filtered_seqs[key] = bs
            
        # TODO: move all print statements to a log file
        # print("Number of sites mutated:", sites_mutated)
        # print("Number of failed mutations:", failed_mutations)
        # print("Number of sites removed:", sites_removed)
        print("Number of sites after filtering:", len(filtered_seqs))
        
        print("Creating m1,m2,m3 sequences...")
        # for each of the filtered sequence, create m1,m2,m3 sequences
        seqdict = {}
        funcdict = {}
        for key in filtered_seqs:
            # Visualization part
            seqdict["%s-wt" % key] = filtered_seqs[key].sequence
            bs = filtered_seqs[key]
            # get m1,m2,m3 for each wt
            for idx, mut in enumerate([[0], [1], [0, 1]]):
                # here we mutate on the first, second, and both sites
                # mut is the index of the site to abolish
                to_remove = bs.remove_pos(mut)
                mutseq = bs.abolish_sites(to_remove, mode="to_eliminate", 
                                          escore_threshold=mutate_cutoff)
                seqdict["%s-m%d" % (key, idx + 1)] = mutseq.sequence
                funcdict["%s-m%d" % (key, idx + 1)] = mutseq.plot_functions
         
            for e in list(itertools.product(egaps, thresholds)): 
                egapthres = e[0]
                ecutoff = e[1]
                # check that wt, m1, m2, m3 are valid
                if coopfilter.check_all_seqs(seqdict["%s-wt"%key], seqdict["%s-m1"%key],
                                             seqdict["%s-m2"%key], seqdict["%s-m3"%key],
                                             filtered_seqs[key].get_sites_dict(), escores,
                                             escore_cutoff=ecutoff, escore_gap=egapthres):
                    bsites_dict = filtered_seqs[key].get_sites_dict()
                    filtered_probes.append({"key":key,
                                            "wt": seqdict["%s-wt"%key],
                                            "m1": seqdict["%s-m1"%key],
                                            "m2": seqdict["%s-m2"%key],
                                            "m3": seqdict["%s-m3"%key],
                                            "core1_start": bsites_dict["core_start_1"],
                                            #"core1_mid": bsites_dict["core_mid_1"],
                                            "core1_end": bsites_dict["core_end_1"],
                                            "site1_pref": bsites_dict["score_1"],
                                            "core2_start": bsites_dict["core_start_2"],
                                            #"core2_mid": bsites_dict["core_mid_2"],
                                            "core2_end": bsites_dict["core_end_2"],
                                            "site2_pref": bsites_dict["score_2"],
                                            "ecutoff": ecutoff,
                                            "egapthres": egapthres,
                                            "distance": filtered_seqs[key].get_sites_dist(),
                                            "sites_in_peak": sitenum,
                                            "peak_length": peaklen
                                            })
                    break # the sequence passes the filtering check, so stop
        filtered_es_preds = {}
        filtered_esplots = {}
        filtered_model_preds = {}
        filtered_model_plots = {}
        for protein in proteins:
            protein_num = proteins.index(protein)
            filtered_es_preds[protein] = escores[protein].predict_sequences(seqdict, key_colname="key")
            filtered_esplots[protein] = escores[protein].make_plot_data(filtered_es_preds[protein], color=colors[protein_num][0])

            filtered_model_preds[protein] = models[protein].predict_sequences(seqdict, 
                                                                     key_colname="key", 
                                                                     predict_flanks=predict_flanks)
            filtered_model_plots[protein] = models[protein].make_plot_data(filtered_model_preds[protein], 
                                                                  color=colors[protein_num][1],
                                                                  show_model_flanks=show_model_flanks) 
        sp.plot_seq_combine([filtered_esplots, filtered_model_plots], 
                            filepath="%splot_mut_d%d_p%d.pdf" % (analysis_path,sitenum,peaklen))

        break
    # initialize a dataframe of filtered probes with m1, m2, m3
    fp_df = pd.DataFrame(filtered_probes)
    req_cols = ["key", "wt", "m1", "m2", "m3", "flank_left", "flank_right", 
                "core1_start", "core1_mid", "core1_end", "core2_start", "core2_mid",
                "core2_end", "ecutoff", "egapthres", "distance", "sites_in_peak",
                "peak_length", "coordinate"]

    df_init = seqdf[["key","flank_left","flank_right","chr","seqstart","seqend"]] # will produce warning but it is fine
    df_init["coordinate"] = df_init.apply(lambda row: "%s:%d-%d" % (row["chr"],row["seqstart"],row["seqend"]), axis=1)
    df_init.drop(['chr', 'seqstart', 'seqend'], axis=1)
    df_out = pd.concat([fp_df.set_index('key'),df_init.set_index('key')], axis=1, join='inner').reset_index()

    if fp_df.columns.isin(req_cols).all(): # this condition check if the columns match exactly
        print("Saving result in %s/mutated_probes.tsv" % analysis_path)
        df_out.to_csv("%s/mutated_probes.tsv" % (analysis_path),sep="\t",index=False,columns=req_cols)
    else:
        print("No mutation rows found in %s" % sitepath)
