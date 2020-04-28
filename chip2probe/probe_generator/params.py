paramlist = {
    # path to all_sites.tsv
    "analysis_path": "",
    # provide either escore_long_path OR (escore_short_path AND escore_map_path)
    # Download escore short and map from Qbic http://qbic.genome.duke.edu/downloads
    # next to 'Download all e-score tables used in our prediction'
    # read escore file
    "escore_short_paths": {'ets1': "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/ets1_escores.txt"},
    "escore_map_paths": {'ets1': "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/index_short_to_long.csv"},

    # path will be provided by user
    "escore_long_paths": {'ets1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/escores/Ets1_8mers_11111111.txt",
                          'runx1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/escores/Runx1_8mers_11111111.txt"
                          },
    # set the model used, either kompas or imads
    "sitecall": "imads",
    "kmer_align_paths": {'ets1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/kompas/ets1/Ets1_kmer_alignment.txt",
                         'runx1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/kompas/runx1/Runx1_kmer_alignment.txt"
                         },
    "modelcores": {'ets1': ["GGAA", "GGAT"],
                   'runx1': ["GAGGT", "GCGGC", "GCGGG", "GCGGT", "GTGGC", "GTGGG", "GTGGT"]
                   },
    "modelpaths": {'ets1': ["/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
                            "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"],
                   'runx1': ["/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GAGGT_1a2a3mer_format.model",
                             "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGC_1a2a3mer_format.model",
                             "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGG_1a2a3mer_format.model",
                             "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGT_1a2a3mer_format.model",
                             "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGC_1a2a3mer_format.model",
                             "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGG_1a2a3mer_format.model",
                             "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGT_1a2a3mer_format.model"]
                   },
    "imads_bed_path": "/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_files/predictions/hg19_0005_Ets1_filtered.bed",
    "imads_cutoff": {'ets1': 0.2128,
                     'runx1': 0.3061},
    "model_kmers" : [1,2,3],
    "imads_sitelen": 20,
    "proteins": ['ets1'],
    # list of escore gaps
    "egaps": [0],
    # list of escore thresholds
    "thresholds": [0.4],
    # cutoff for significant escore. If mutate is True, any non-centered sites below this cutoff are
    # mutated and removed
    "mutate_cutoff": 0.38,
    # gap used in determining what to mutate
    "mutate_gap": 0,
    # imads would usually use flanks to predict
    "predict_flanks": True,
    # boolean whether to show flanks used by the model. Kompas has no flanks, only imads does
    "show_model_flanks": False,
    # user can determine colors for the proteins as a list of tuples (escore color, model color)
    "colors": [],
    "generate_plots": False,
    # indicate whether or not we require all m1,m2,m3 sequences for each wild type
    "get_complete_mutated": True,
    # 1: make weak sites, only for imads
    # 2: bring b1 closer to b2, fixing b2, adding flank
    # 3: bring b2 closer to b1, fixing b1, adding flank
    # 4: bring b1 closer to b2, fixing b2, not adding flank
    # 5: fix distance, move binding sites to the right, adding flank
    # 6: fix distance, move binding sites to the left, adding flank
    # 7: fix distance, move binding sites to the right, not adding flank
    "customs": [2, 3, 4, 5, 6, 7],
    "tf_source": "ets1_k562",
    # number of custom sequences
    "ncustom": 750,
    # primer sequence to be used
    "primer": "GTCTTGATTCGCTTGACGCTGCTG",
    # max number of nucelotides to be mutated at the junctions for cleaning
    "max_mutate_count": 2
}
