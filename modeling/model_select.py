import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL
sys.path.append("/Users/faricazjj/Desktop/homotf/chip2probe")
import os
sys.path.append(os.path.join('C:/', 'Users', 'Farica Zhuang', 'Desktop', 'homotf', 'chip2probe'))
from trainingdata.training import Training
import util.util as util
from best_model import BestModel
from trainingdata.dnashape import DNAShape

from sklearn import ensemble, model_selection, metrics, tree
import pandas as pd
import numpy as np
import scipy
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import pickle

import subprocess

def get_numeric_label(training):
    # hard coded but change add to anti coop / additive when needed
    train = training['label'].map({'cooperative': 1, 'additive': 0})
    return train


def display_output(fpr_list, tpr_dict, auc_dict, path, title):
    """
        This plots the average ROC curve of all the classifiers in a single plot
    """
    plt.clf() # first, clear the canvas
    plt.plot([0, 1], [0, 1], color="gray", alpha=0.5, lw=0.3)#linestyle="--",
    for key in  tpr_dict:
        tpr_list = tpr_dict[key]
        auc = auc_dict[key]

        #for key in fpr_dict:
        """
        ln = key.split(",")
        if len(ln) == 3:
            plt.plot(fpr_list, tpr_list, lw=2, label='%s: %.3f' % (key,auc))
        elif len(ln) == 2:
            plt.plot(fpr_list, tpr_list, lw=1, label='%s: %.3f' % (key,auc))
        else:
        plt.plot(fpr_list, tpr_list, linestyle="--", lw=1, label='%s: %.3f' % (key,auc))
        """
        plt.plot(fpr_list, tpr_list, lw=2, label='%s: %.2f' % (key,auc))

        # Show the ROC curves for all classifiers on the same plot
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')

    plt.title(title)
    leg = plt.legend(loc="lower right",title="AUC for each combination:")
    leg._legend_box.align = "left"
    plt.savefig(path, dpi=1200)

def plot_auc(x_train_dict, title="Average ROC Curves", plotname="auc.png"):
    tpr_dict = {key:[] for key in x_train_dict}
    auc_dict = {key:[] for key in x_train_dict}
    acc_dict = {key:[] for key in x_train_dict}

    fpr_lim=100
    base_fpr = np.linspace(0, 1, 101)[:fpr_lim+1]
    cv = model_selection.KFold(n_splits=10,shuffle=True)

    random_x = next(iter(x_train_dict.values()))[0]
    x_train = random_x.values.tolist() # .loc[:,random_x.columns != 'label']
    y_train = random_x['label'].values

    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=10,random_state=0)
    for train_idx,test_idx in cv.split(x_train):
        tprs = []
        aucs_val = []
        # compare with just distance
        for key in x_train_dict:
            # need to convert this with index, somehow cannot do
            # x_train[train_idx] for multi features
            xt_df = x_train_dict[key][0]
            xt = xt_df.loc[:,xt_df.columns != 'label'].values.tolist() # need to be a list of list

            data_train = [xt[i] for i in train_idx]
            data_test = [xt[i] for i in test_idx]
            lbl_train = [y_train[i] for i in train_idx]
            lbl_test = [y_train[i] for i in test_idx]

            model = x_train_dict[key][1]

            model = model.fit(data_train, lbl_train)
            y_score = model.predict_proba(data_test)
            fpr, tpr, _ = metrics.roc_curve(lbl_test, y_score[:, 1])
            #auc = metrics.roc_auc_score(lbl_test, y_score[:,1])
            tpr = scipy.interp(base_fpr, fpr, tpr) # add points to the plotting
            res_auc = metrics.auc(base_fpr, tpr)
            tpr_dict[key].append(np.insert(tpr,0,0)) # insert 0 so we start the plot from the bottom left
            auc_dict[key].append(res_auc)

            y_pred = model.predict(data_test)
            acc_dict[key].append(accuracy_score(lbl_test, y_pred))

    mean_tpr = {k:np.array(tpr_dict[k]).mean(axis=0) for k in tpr_dict}
    mean_auc = {k:np.array(auc_dict[k]).mean(axis=0) for k in auc_dict}
    mean_acc = {k:np.array(acc_dict[k]).mean(axis=0) for k in acc_dict}
    print("Mean accuracy", mean_acc, "\nMean auc", mean_auc)

    # left append 0 in base fpr just so we start at 0 (we did the same for the tpr)
    base_fpr = np.insert(base_fpr,0,0)
    display_output(base_fpr, mean_tpr, mean_auc, path=plotname, title=title)

if __name__ == "__main__":
    trainingpath = "train1.tsv"
    trainingpath = "trainingdata/training_new.csv"
    #shapepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/dnashape/training_p01_adjusted_reversed"
    shapepath = "training_p01_adjusted_reversed"
    
    ds = DNAShape(shapepath)

    rf_param_dict = {
                    'n_estimators': [200,500,1000,1500],
    				'max_depth': [5,10,15],
                    'min_samples_split': [10,20,30]
    			}
    dt_param_dict = {
    				"criterion" : ['gini', 'entropy'],
                    "min_samples_split" : [i for i in range(2,41)],
                    "min_samples_leaf" : [i for i in range(2,31)]
    			}

    df = pd.read_csv(trainingpath, sep=",")
    #df = df[df["orientation"] == "TT"]

    t = Training(df, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])

    # xtr = {
    #         "all":
    #             BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               "distance":{"type":"numerical"},
    #                               "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True},
    #                               "sitepref": {},
    #                               "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
    #                                "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
    #                                "flankseq": {"k":3, "seqin":4, "smode":"strength"},
    #                                "flankseq": {"k":3, "seqin":-4, "smode":"strength"}
    #                           })
    #             ).run_all(),

    #         "top10":
    #             BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               "distance":{"type":"numerical"},
    #                               "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True},
    #                               "sitepref": {},
    #                               "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
    #                                "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
    #                                "flankseq": {"k":3, "seqin":4, "smode":"strength"},
    #                                "flankseq": {"k":3, "seqin":-4, "smode":"strength"}
    #                           }),
    #                       topn=10
    #             ).run_all()
    # }
   
    xtr = {
            "distance":
                BestModel(clf="RF",
                          param_dict=rf_param_dict,
                          train_data=t.get_training_df({
                                  "distance":{"type":"numerical"},
                                  #"orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True, "include":"F"}
                              })
                ).run_all(),
            "flank seq":
                BestModel(clf="RF",
                          param_dict=rf_param_dict,
                          train_data=t.get_training_df({
                                  "flankseq": {"k":3, "seqin":4, "smode":"strength"},
                                  "flankseq": {"k":3, "seqin":-4, "smode":"strength"},
                                  # "flankshape": {"ds":ds, "seqin":4, "smode":"positional"},
                                  # "flankshape": {"ds":ds, "seqin":-4, "smode":"positional"},
                                  #"orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True, "include":"F"}
                              })
                ).run_all(),
            "distance, flank seq":
                BestModel(clf="RF",
                          param_dict=rf_param_dict,
                          train_data=t.get_training_df({
                                  "distance":{"type":"numerical"},
                                  "flankseq": {"k":3, "seqin":4, "smode":"strength"},
                                  "flankseq": {"k":3, "seqin":-4, "smode":"strength"},
                                   # "flankshape": {"ds":ds, "seqin":4, "smode":"positional"},
                                   # "flankshape": {"ds":ds, "seqin":-4, "smode":"positional"},
                                   #"orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True, "include":"F"}
                              })
                ).run_all(),
             "top10":
             	BestModel(clf="RF",
                          param_dict=rf_param_dict,
                          train_data=t.get_training_df({
                                  "distance":{"type":"numerical"},
                                  "flankseq": {"k":3, "seqin":4, "smode":"strength"},
                                  "flankseq": {"k":3, "seqin":-4, "smode":"strength"},
                                  # "flankshape": {"ds":ds, "seqin":4, "smode":"positional"},
                                  # "flankshape": {"ds":ds, "seqin":-4, "smode":"positional"},
                                  #"orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True, "include":"F"}
                              }),
                           topn=10
                ).run_all()
        }

   

    # xtr = {
    #         "distance":
    #             BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               "distance":{"type":"numerical"},
    #                               #"orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True, "include" : "F"}
    #                           })
    #             ).run_all(),
    #         # "orientation":
    #         #     BestModel(clf="RF",
    #         #               param_dict=rf_param_dict,
    #         #               train_data=t.get_training_df({
    #         #                       #"distance":{"type":"numerical"},
    #         #                       "orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True}
    #         #                   })
    #         #     ).run_all(),
    #         "strength":
    #             BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               #"distance":{"type":"numerical"},
    #                               #"orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True}
    #                               "sitepref": {}
    #                           })
    #             ).run_all(),
    #         "distance, strength":
    #             BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               "distance":{"type":"numerical"},
    #                               #"orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True}
    #                               "sitepref": {}
    #                           })
    #             ).run_all()
            # "distance, orientation":
            #     BestModel(clf="RF",
            #               param_dict=rf_param_dict,
            #               train_data=t.get_training_df({
            #                       "distance":{"type":"numerical"},
            #                       "orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True}
            #                       #"sitepref": {}
            #                   })
            #     ).run_all(),
            # "distance, strength, orientation":
            #     BestModel(clf="RF",
            #               param_dict=rf_param_dict,
            #               train_data=t.get_training_df({
            #                       "distance":{"type":"numerical"},
            #                       "orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True},
            #                       "sitepref": {}
            #                   })
            #     ).run_all()
       # }
            # "dist-flank-seq":
            #     BestModel(clf="RF",
            #               param_dict=rf_param_dict,
            #               train_data=t.get_training_df({
            #                       "distance":{"type":"numerical"},
            #                       "flankseq": {"k":3, "seqin":4, "smode":"strength"},
            #                       "flankseq": {"k":3, "seqin":-4, "smode":"strength"},
            #                       "orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True}
            #                   }),
            #                       # 
            #                   # })
            #     ).run_all(),
            # "all":
            #     BestModel(clf="RF",
            #               param_dict=rf_param_dict,
            #               train_data=t.get_training_df({
            #                       "distance":{"type":"numerical"},
            #             "sitepref": {},
            #             "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
            #                   })
            #     ).run_all(),
            # "top10-all":
            #     BestModel(clf="RF",
            #               param_dict=rf_param_dict,
            #               train_data=t.get_training_df({
            #                       "distance":{"type":"numerical"},
            #             "sitepref": {},
            #             "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
            #                   },
            #                ),
            #                topn=10
            #     ).run_all()
           

    plot_auc(xtr, "Average ROC Curves Using RF", "dist_flankseq_strength_rf.png")
    """
    # save the first model
    xt = t.get_feature_all({
        "distance":{"type":"numerical"},
        "flankshape": {"ds":ds, "seqin":5, "smode":"positional"},
        "flankshape": {"ds":ds, "seqin":-3, "smode":"positional"}
        })
    topf = ['dist_numeric', 'prot_outer_wk_pos_2', 'helt_outer_str_pos_2', 'prot_outer_str_pos_1',  'roll_outer_wk_pos_2', 'roll_outer_str_pos_1', 'mgw_outer_wk_pos_0', 'prot_outer_str_pos_0', 'prot_outer_str_pos_2', 'mgw_outer_str_pos_0']
    xt_df = pd.DataFrame(xt)#[topf]
    xtlist = xt_df.values.tolist()
    y_train = get_numeric_label(t.df).values
    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=10,random_state=0,min_samples_leaf=20)
    m = rf.fit(xtlist, y_train)
    tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
            feature_names = xt_df.columns,
            class_names = ['additive','cooperative'],
            rounded = True, proportion = False,
            precision = 2, filled = True)
    subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])

    #plot_auc(xtr, y_train, "auc.png")
    # xt = xt_df.values.tolist()
    # rf.fit(xt,y_train)
    # pickle.dump(rf, open("model1.sav", 'wb'))
    """
