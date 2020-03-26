import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL
#sys.path.append("/Users/faricazjj/Desktop/homotf/chip2probe")
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
        plt.plot(fpr_list, tpr_list, lw=2, label='%s: %.3f' % (key,auc))

        # Show the ROC curves for all classifiers on the same plot
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')

    plt.grid()
    plt.title(title)
    leg = plt.legend(loc="lower right",title="AUC for each combination:")
    leg._legend_box.align = "left"
    plt.savefig(path)

def plot_auc(x_train_dict, y_train, title="Average ROC Curves", plotname="auc.png"):
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

            model = rf# x_train_dict[key][1]

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
    #trainingpath = "trainingdata/training_new.csv"
    shapepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/dnashape/training_p01_adjusted_reversed"
    ds = DNAShape(shapepath)

    rf_param_dict = {
                    'n_estimators': [i for i in range(100,601,100)],
    				'max_depth': [i for i in range(3,11)]
    			}
    dt_param_dict = {
    				"criterion" : ['gini', 'entropy'],
                    "min_samples_split" : [i for i in range(20,41)],
                    "min_samples_leaf" : [i for i in range(20,31)]
    			}
    df = pd.read_csv(trainingpath, sep=",")
    #df = df[dftrain["orientation"] == "HT/TH"]

    t = Training(df, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])
    y_train = get_numeric_label(t.df).values

    xtr = {
            "dist-ori":
                BestModel(clf="RF",
                          param_dict=rf_param_dict,
                          train_data=t.get_training_df({
                                  "distance":{"type":"numerical"},
                                  "orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True}
                              })
                ).run_all(),
            "ori-flank-seq":
                BestModel(clf="RF",
                          param_dict=rf_param_dict,
                          train_data=t.get_training_df({
                                  "flankseq": {"k":3, "seqin":4, "smode":"strength"},
                                  "flankseq": {"k":3, "seqin":-4, "smode":"strength"},
                                  "orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True}
                              }),
                                  # 
                              # })
                ).run_all(),
            "dist-ori-flank-seq":
            	BestModel(clf="RF",
                          param_dict=rf_param_dict,
                          train_data=t.get_training_df({
            					  "distance":{"type":"numerical"},
                                  "flankshape": {"ds":ds, "seqin":4, "smode":"strength"},
                                  "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
                                  "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                              })
                ).run_all(),
            # "flankshape":
            #     BestModel(clf="RF",
            #               param_dict=rf_param_dict,
            #               train_data=t.get_training_df({
            #                       "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
            #                       "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
            #                       "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
            #                   })
            #     ).run_all(),
            # "dist-ori-flankshape":
            #     BestModel(clf="RF",
            #               param_dict=rf_param_dict,
            #               train_data=t.get_training_df({
            #                       "distance":{"type":"numerical"},
                                   # "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
                                   # "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
                                   # "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                #               }),
                # ).run_all(),
             "top10":
             	BestModel(clf="RF",
                          param_dict=rf_param_dict,
                          train_data=t.get_training_df({
                                  "distance":{"type":"numerical"},
                                  "flankseq": {"k":3, "seqin":4, "smode":"strength"},
                                  "flankseq": {"k":3, "seqin":-4, "smode":"strength"},
                                  "orientation":{"positive_cores":["GGAA", "GGAT"], "one_hot":True}
                                  # "flankshape": {"ds":ds, "seqin":4, "smode":"positional"},
                                  # "flankshape": {"ds":ds, "seqin":-3, "smode":"positional"}
                                  # "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
                                  # "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
                                  # "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
                              }),
                           topn=10
                ).run_all()
        }


    plot_auc(xtr, y_train, "Average ROC Curves Using RF for All Orientations", "dist_flank_seq_auc.png")

    # # save the first model
    # dt = tree.DecisionTreeClassifier(min_samples_split=27, min_samples_leaf=25, criterion="entropy")
    # #plot_auc(xtr, y_train, "auc.png")
    # xt = t.get_feature_all({
    #     "distance":{"type":"numerical"}
    #     })
    # xt = pd.DataFrame(xt).values.tolist()
    # dt.fit(xt,y_train)
    # pickle.dump(rf, open("model1.sav", 'wb'))
