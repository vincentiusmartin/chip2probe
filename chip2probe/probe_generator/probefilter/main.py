'''
Created on Sep 30, 2019

@author: vincentiusmartin
'''

# these paths need to be handled better
import sys
import pandas as pd
import itertools
import os
from chip2probe.probe_generator.probefilter import sitespredict
from sitespredict import *
from cooperative import coopfilter
from chip2probe.util import bio as bio

if __name__ == '__main__':
    # path to all_sites.tsv
    #analysis_path = sys.argv[2]
    # provide either escore_long_path OR (escore_short_path AND escore_map_path)
    # Download escore short and map from Qbic http://qbic.genome.duke.edu/downloads
    # next to 'Download all e-score tables used in our prediction'
    # read escore file
    escore_short_paths = {'ets1': "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/ets1_escores.txt"}
    escore_map_paths = {'ets1': "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/index_short_to_long.csv"}

    escore_long_paths = {'ets1': "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/Ets1_8mers_11111111.txt"}

    #curdir = "/Users/vincentiusmartin/Research/chip2gcPBM/"
    #analysis_path = curdir + "result/cistrome_ets1_37927/analysis_result/"

    # set the model used, either kompas or imads
    method = "imads"
    proteins = list(escore_long_paths.keys())
    print(proteins)

    # # create escore object for each protein
    # escores = {}
    # escores['ets1'] = PBMEscore(escore_long_path)

    # models = {}
    # # initialize binding site prediction model used based on method
    # if method=='imads':
    #     # read imads prediction
    #     modelcores = ["GGAA", "GGAT"]
    #     modelpaths = ["/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
    #                   "/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"]
    #     imads_models = [iMADSModel(modelpath, modelcore, 20, [1, 2, 3])
    #                 for modelpath, modelcore in zip(modelpaths, modelcores)]
    #     models['ets1'] = iMADS(imads_models, 0.2128)  # 0.2128 is for the ETS1 cutoff
    # elif method == 'kompas':
    #     models['ets1'] = Kompas()
    # else:
    #     raise Exception("Supported methods are imads and kompas. Got {}".format(method))

    # pd.set_option('display.max_columns', None)
    # sitepath = "%s/sites_all.tsv" % analysis_path
    # seqdf = pd.read_csv(sitepath, sep='\t')

    # spcomb = [(x, y) for x in seqdf["sites.in.peak"].unique() for y in seqdf["peaklen"].unique()]
    # filtered_probes = []  # we will fill this
    # for comb in spcomb:
    #     sitenum = comb[0]
    #     peaklen = comb[1]
    #     print("Working on sitenum %d, peaklen %d" % (sitenum,peaklen))
    #     curdf = seqdf.loc[(seqdf["sites.in.peak"] == sitenum) & (seqdf["peaklen"] == peaklen)]
    #     for protein in proteins:
    #         es_preds = escores[protein].predict_sequences(curdf, key_colname="key")
    #         esplots = escores[protein].make_plot_data(es_preds)

    #         model_preds = imads.predict_sequences(curdf, key_colname="key")
    #         imads_plots = imads.make_plot_data(imads_preds)

    #     sp = SitesPlotter()
    #     # if need to plot, uncomment this
    #     sp.plot_seq_combine([esplots,imads_plots], filepath="%s/sitesplot_d%d_p%d.pdf" % (analysis_path,sitenum,peaklen))

    #     # get filtered sequences
    #     filtered_seqs = {}
    #     flanks = {}
    #     print("Site filtering...")
    #     # get sequences with 2 significant binding sites
    #     for key in es_preds:
    #         bs = Sequence(es_preds[key], model_preds[key], escore_cutoff=0.35)
    #         if bs.site_count() == 2:
    #             filtered_seqs[key] = bs

    #     seqdict = {}
    #     funcdict = {}
    #     for key in filtered_seqs:
    #     #for key in ["sequence11"]:
    #         # Visualization part
    #         seqdict["%s-wt" % key] = filtered_seqs[key].sequence
    #         for idx,mut in enumerate([[0],[1],[0,1]]):
    #             # here we mutate on the first, second, and both sites
    #             # mut is the index of the site to abolish
    #             mutseq = filtered_seqs[key].abolish_sites(mut,escore)
    #             seqdict["%s-m%d" % (key,idx + 1)] = mutseq.sequence
    #             funcdict["%s-m%d" % (key,idx + 1)] = mutseq.plot_functions
    #         for e in list(itertools.product([0,1,2],[0.41, 0.4, 0.35])): # ,0.35
    #             egapthres = e[0]
    #             ecutoff = e[1]
    #             if coopfilter.filter_coopseq(seqdict["%s-wt"%key], seqdict["%s-m1"%key],
    #                                  seqdict["%s-m2"%key], seqdict["%s-m3"%key],
    #                                  filtered_seqs[key].get_sites_dict(), escore,
    #                                  escore_cutoff=ecutoff, escore_gap = egapthres):
    #                 bsites_dict = filtered_seqs[key].get_sites_dict()
    #                 filtered_probes.append({"key":key,
    #                                     "wt":seqdict["%s-wt"%key],
    #                                     "m1":seqdict["%s-m1"%key],
    #                                     "m2":seqdict["%s-m2"%key],
    #                                     "m3":seqdict["%s-m3"%key],
    #                                     "core1_start":bsites_dict["core_start_1"],
    #                                     "core1_end":bsites_dict["core_end_1"],
    #                                     "site1_pref":bsites_dict["imads_score_1"],
    #                                     "core2_start":bsites_dict["core_start_2"],
    #                                     "core2_end":bsites_dict["core_end_2"],
    #                                     "site2_pref":bsites_dict["imads_score_2"],
    #                                     "ecutoff": ecutoff,
    #                                     "egapthres":egapthres,
    #                                     "distance":filtered_seqs[key].get_sites_dist(),
    #                                     "sites_in_peak":sitenum,
    #                                     "peak_length":peaklen
    #                                     })
    #                 break # the sequence passes the filtering check, so stop
    #     pp = escore.make_plot_data(escore.predict_sequences(seqdict),additional_functions=funcdict)
    #     sp.plot_seq_combine([pp], filepath="%s/plot_mut_d%d_p%d.pdf" % (analysis_path,sitenum,peaklen))

    # # initialize a dataframe of filtered probes with m1, m2, m3
    # fp_df = pd.DataFrame(filtered_probes)
    # req_cols = ["key", "wt", "m1", "m2", "m3", "flank_left", "flank_right", 
    #             "core1_start", "core1_mid", "core1_end", "core2_start", "core2_mid",
    #             "core2_end", "ecutoff", "egapthres", "distance", "sites_in_peak",
    #             "peak_length", "coordinate"]

    # df_init = seqdf[["key","flank_left","flank_right","chr","seqstart","seqend"]] # will produce warning but it is fine
    # df_init["coordinate"] = df_init.apply(lambda row: "%s:%d-%d" % (row["chr"],row["seqstart"],row["seqend"]), axis=1)
    # df_init.drop(['chr', 'seqstart', 'seqend'], axis=1)
    # df_out = pd.concat([fp_df.set_index('key'),df_init.set_index('key')], axis=1, join='inner').reset_index()

    # if fp_df.columns.isin(req_cols).all(): # this condition check if the columns match exactly
    #     print("Saving result in %s/mutated_probes.tsv" % analysis_path)
    #     df_out.to_csv("%s/mutated_probes.tsv" % (analysis_path),sep="\t",index=False,columns=req_cols)
    # else:
    #     print("No mutation rows found in %s" % sitepath)
