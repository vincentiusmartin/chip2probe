import os
os.chdir("../..")

import pandas as pd

from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.bestmodel import BestModel
import chip2probe.modeler.plotlib as pl
from sklearn import ensemble, tree
import subprocess

# "shape":
#     BestModel(clf="sklearn.ensemble.RandomForestClassifier",
#       param_grid=rf_param_grid,
#       train_data=ct.get_training_df({
#             "distance":{"type":"numerical"},
#             "shape_in":{"seqin":3, "poscols":['ets_pos','runx_pos']},
#             "shape_out":{"seqin":-2, "poscols":['ets_pos','runx_pos']}
#         })
#     ).run_all()

if __name__ == "__main__":
    trainingpath = "output/heterotypic/EtsRunx_v1/ch1_ch2/training_pwm.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    df['label'] = df['label'].replace('independent', 'additive')
    ct = CoopTrain(df)
    pd.set_option("display.max_columns",None)

    rf_param_grid = {
        'n_estimators': [500], #[500,750,1000],
        'max_depth': [10], #[5,10,15],
        "min_samples_leaf": [10], #[5,10,15],
        "min_samples_split" : [20], #[5,10,15]
    }

    best_models = {
        "distance,orientation":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}
                })
            ).run_all(),
        "distance,orientation,sequence":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}},
                    "sequence_in":{"seqin":5, "poscols":['ets_pos','runx_pos'], "namecol":"Name"},
                    "sequence_out":{"seqin":-3, "poscols":['ets_pos','runx_pos'], "namecol":"Name"}
                })
            ).run_all(),
        "distance,orientation,shape":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}},
                    "shape_in":{"seqin":5, "poscols":['ets_pos','runx_pos']},
                    "shape_out":{"seqin":-2, "poscols":['ets_pos','runx_pos']}
                })
            ).run_all(),
        "distance,orientation,sequence,shape":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}},
                    "shape_in":{"seqin":5, "poscols":['ets_pos','runx_pos']},
                    "shape_out":{"seqin":-2, "poscols":['ets_pos','runx_pos']},
                    "sequence_in":{"seqin":5, "poscols":['ets_pos','runx_pos'], "namecol":"Name"},
                    "sequence_out":{"seqin":-3, "poscols":['ets_pos','runx_pos'], "namecol":"Name"}
                })
            ).run_all(),
        "distance,orientation,strength":
            BestModel(clf="sklearn.ensemble.RandomForestClassifier",
              param_grid=rf_param_grid,
              train_data=ct.get_training_df({
                    "distance":{"type":"numerical"},
                    "affinity": {"colnames": ("ets_score","runx_score")},
                    "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}
                })
            ).run_all(),
    }

    pl.plot_model_metrics(best_models, cvfold=10, score_type="auc", varyline=True, title="Average ROC Curves for Ets1-Runx1\n(using shape and sequence features)")

    # feature_dict = {
    #     "distance":{"type":"numerical"},
    #     "affinity": {"colnames": ("ets_score","runx_score")},
    #     "orientation": {"relative":False, "pos_cols": {"ets_pos":"ets_ori", "runx_pos":"runx_ori"}}
    # }
    # train = ct.get_feature_all(feature_dict)
    # label = ct.get_numeric_label({'cooperative': 1, 'additive': 0})
    # rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=5, min_samples_leaf=10, min_samples_split=10)
    # m = rf.fit(train.values.tolist(),label)
    #
    # tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
    #         feature_names = train.columns,
    #         class_names = ['additive','cooperative'],
    #         rounded = True, proportion = False,
    #         precision = 2, filled = True)
    # subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])
