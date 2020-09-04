import pandas as pd

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pbmescore import PBMEscore
from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.sitesplotter import SitesPlotter

seq = "TTACGGCAAGCGGGCCGGAAGCCACTCCTCGAGTCT"
df = pd.DataFrame([[seq, 'seq1']], columns=['sequence', 'key'])

mutate_cutoff = 0.38

# path will be provided by user
escore_long_paths= {'ets1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/escores/Ets1_8mers_11111111.txt",
                      'runx1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/escores/Runx1_8mers_11111111.txt"
                      }
kmer_align_paths= {'ets1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/kompas/ets1/Ets1_kmer_alignment.txt",
                     'runx1': "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/kompas/runx1/Runx1_kmer_alignment.txt"
                     }
modelcores= {'ets1': ["GGAA", "GGAT"],
               'runx1': ["GAGGT", "GCGGC", "GCGGG", "GCGGT", "GTGGC", "GTGGG", "GTGGT"]
               }
modelpaths= {'ets1': ["/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
                        "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"],
               'runx1': ["/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GAGGT_1a2a3mer_format.model",
                         "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGC_1a2a3mer_format.model",
                         "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGG_1a2a3mer_format.model",
                         "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GCGGT_1a2a3mer_format.model",
                         "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGC_1a2a3mer_format.model",
                         "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGG_1a2a3mer_format.model",
                         "/Users/faricazjj/Box/homotf/chip2probe/chip2probe/data/imads/runx1/Runx1_10nM_Bound_filtered_normalized_logistic_transformed_20bp_GTGGT_1a2a3mer_format.model"]
               }
imads_cutoff= {'ets1': 0.2128,
                 'runx1': 0.3061}
escores = {}
models = {}
proteins = ['ets1', 'runx1']

for tf in proteins:
    escores[tf] = PBMEscore(escore_long_paths[tf])

    models[tf] = Kompas(protein=tf, threshold=mutate_cutoff,
                                 kmer_align_path=kmer_align_paths[tf])

es_preds = {}
esplots = {}
model_preds = {}
model_plots = {}
colors = [('crimson', 'plum'), ('steelblue', 'lightblue')]

# initialize escore and model objects for each protein
for protein in proteins:
    protein_num = proteins.index(protein)
    es_preds[protein] = escores[protein].predict_sequences(df, key_colname="key")
    esplots[protein] = escores[protein].make_plot_data(es_preds[protein], color=colors[protein_num][0])

    model_preds[protein] = models[protein].predict_sequences(df, key_colname="key", predict_flanks=False)
    model_plots[protein] = models[protein].make_plot_data(model_preds[protein],
                                                          color=colors[protein_num][1],
                                                          show_model_flanks=False)

# Generate plots

sp = SitesPlotter()
# if need to plot, uncomment this
sp.plot_seq_combine([esplots, model_plots],
                    filepath="outlier.pdf")
