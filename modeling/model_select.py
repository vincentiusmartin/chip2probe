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
import matplotlib.pyplot as plt
import pickle

import subprocess

def get_numeric_label(training):
    # hard coded but change add to anti coop / additive when needed
    train = training['label'].map({'cooperative': 1, 'additive': 0})
    return train


def display_output(xy, score_dict, path, title, score_type="auc"):
    """
        This plots the average ROC curve of all the classifiers in a single plot
    """
    plt.clf() # first, clear the canvas
    # if score_type == "pr":
    #     plt.plot([0, 1], [1, 0], color="gray", alpha=0.5, lw=0.3)#linestyle="--",
    # else:
    if score_type == "auc":
        plt.plot([0, 1], [0, 1], color="gray", alpha=0.5, lw=0.3)#linestyle="--",
    for key in  xy:
        score = score_dict[key]
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
        plt.plot(xy[key]['x'], xy[key]['y'], lw=2, label='%s: %.3f' % (key,score))

        # Show the ROC curves for all classifiers on the same plot
        if score_type == "pr":
            plt.xlabel('Recall')
            plt.ylabel('Precision')
        else:
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')

    plt.grid()
    plt.title(title)
    leg = plt.legend(loc="lower right",title="%s for each combination:"%str.upper(score_type))
    leg._legend_box.align = "left"
    plt.savefig(path)

# plot_performance
def plot_metrics(x_train_dict, title="Average ROC Curves",
                 plotname="auc.png", n_splits=10, score_type="auc"):
    yprob_dict = {key:[] for key in x_train_dict}
    ytrue_dict = {key:[] for key in x_train_dict}
    ypred_dict = {key:[] for key in x_train_dict}
    # auc_dict = {key:[] for key in x_train_dict}
    # acc_dict = {key:[] for key in x_train_dict}

    x_lim=100
    #base_x = np.linspace(0, 1, 101)
    cv = model_selection.KFold(n_splits=n_splits,shuffle=True)

    random_x = next(iter(x_train_dict.values()))[0]
    xt_all = random_x.values.tolist() # .loc[:,random_x.columns != 'label']
    labels = random_x['label'].values

    for train_idx,test_idx in cv.split(xt_all):
        aucs_val = []
        # compare with just distance
        for key in x_train_dict:
            # need to convert this with index, somehow cannot do
            # x_train[train_idx] for multi features
            xt_df = x_train_dict[key][0]
            xt = np.array(xt_df.loc[:,xt_df.columns != 'label'].values.tolist()) # need to be a list of list

            x_train, x_test = xt[train_idx], xt[test_idx]
            y_train, y_test = labels[train_idx], labels[test_idx]

            model = x_train_dict[key][1]
            model = model.fit(x_train, y_train)
            y_score = model.predict_proba(x_test)
            # if score_type == "auc":
            #     x, y, _ = metrics.roc_curve(lbl_test, y_score[:, 1]) #fpr, tpr, threshold
            # elif score_type == "pr":
            #     y, x, _ = metrics.precision_recall_curve(lbl_test, y_score[:, 1]) # precision, recall, threshold
            #auc = metrics.roc_auc_score(lbl_test, y_score[:,1])
            #y = scipy.interp(base_x, x, y) # add points to the plotting
            # res_auc = metrics.auc(x, y) #base_x
            # auc_dict[key].append(res_auc)

            ytrue_dict[key].append(y_test)
            yprob_dict[key].append(y_score[:, 1])
            ypred_dict[key].append(model.predict(x_test))

            # y_pred = model.predict(data_test)
            # acc_dict[key].append(accuracy_score(lbl_test, y_pred))

    ytrue_dict = {k:np.concatenate(ytrue_dict[k]) for k in ytrue_dict}
    yprob_dict = {k:np.concatenate(yprob_dict[k]) for k in yprob_dict}
    ypred_dict = {k:np.concatenate(ypred_dict[k]) for k in ypred_dict}

    if score_type == "auc":
        mets = {k: metrics.roc_curve(ytrue_dict[k], yprob_dict[k]) for k in ytrue_dict} #fpr, tpr, threshold
        for k in mets:
            mets[k] = {"x":mets[k][0], "y":mets[k][1]} #"fpr, tpr"
    elif score_type == "pr":
        mets = {k: metrics.precision_recall_curve(ytrue_dict[k], yprob_dict[k]) for k in ytrue_dict} # precision, recall, threshold
        for k in mets:
            mets[k] = {"x":mets[k][1], "y":mets[k][0]} #recall, precision

    if score_type == "auc":
        scoreval = {k:metrics.auc(mets[k]['x'],mets[k]['y']) for k in ytrue_dict} #{k:np.array(auc_dict[k]).mean(axis=0) for k in auc_dict}
    else: # average_precision_score
        scoreval = {k:metrics.average_precision_score(ytrue_dict[k],yprob_dict[k]) for k in ytrue_dict}
    acc = {k:metrics.accuracy_score(ytrue_dict[k],ypred_dict[k]) for k in ytrue_dict}
    confmat = {k:metrics.confusion_matrix(ytrue_dict[k],ypred_dict[k]).ravel() for k in ytrue_dict}
    print("Mean accuracy", acc, "\nMean %s" % score_type, scoreval)
    print("Confusion_matrix (tn,fp,fn,tp):", confmat)
    # left append 0 in base fpr just so we start at 0 (we did the same for the tpr)
    # if score_type == "auc":
    #     base_x = np.insert(base_x,0,0)
    #     for k in y_dict:
    #         mean_y[key] = np.insert(mean_y[key],0,0)
    display_output(mets, scoreval, path=plotname, title=title, score_type=score_type)

if __name__ == "__main__":
    trainingpath = "train1.tsv"
    #trainingpath = "trainingdata/training_new.csv"
    shapepath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/dnashape/training_p01_adjusted_reversed"
    ds = DNAShape(shapepath)
    score_type = "auc"

    rf_param_dict = {
                    'n_estimators': [500, 1000, 1500],
    				'max_depth': [5, 10, 15],
                    "min_samples_leaf" : [10,15,20],
                    "min_samples_split" : [10,15,20]
    			}
    dt_param_dict = {
    				"criterion" : ['gini', 'entropy'],
                    "min_samples_split" : [i for i in range(20,41)],
                    "min_samples_leaf" : [i for i in range(20,31)]
    			}

    df = pd.read_csv(trainingpath, sep="\t")
    #df = df[df["orientation"] == "TT"]

    t = Training(df, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])


    # xtr = {
    #         "distance-ori":
    #             BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               "distance":{"type":"numerical"},
    #                               "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    #                           })
    #             ).run_all(score_type=score_type),
    #         "flankshape-ori":
    #             BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
    #                               "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
    #                               "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    #                           })
    #             ).run_all(),
    #         "dist-flankshape-ori":
    #             BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               "distance":{"type":"numerical"},
    #                                "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
    #                                "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
    #                                "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    #                           }),
    #             ).run_all(score_type=score_type),
    #          "top10":
    #          	BestModel(clf="RF",
    #                       param_dict=rf_param_dict,
    #                       train_data=t.get_training_df({
    #                               "distance":{"type":"numerical"},
    #                               "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
    #                               "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
    #                               "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    #                           }),
    #                        topn=10
    #             ).run_all()
    #     }
    #
    # plot_metrics(xtr, "Average ROC Curves Using RF for All Orientations", "dist_flank_seq_auc.png",score_type=score_type)

    # save the first model
    xt = t.get_feature_all({
        "distance":{"type":"numerical"},
        "flankshape": {"ds":ds, "seqin":5, "smode":"strength"},
        "flankshape": {"ds":ds, "seqin":-3, "smode":"strength"},
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True}
    })
    topf = ['dist_numeric', 'roll_outer_wk_pos_2', 'helt_outer_wk_pos_2', 'roll_outer_wk_pos_1', 'helt_outer_wk_pos_0', 'roll_outer_wk_pos_0', 'roll_outer_str_pos_0', 'helt_outer_wk_pos_1', 'prot_outer_str_pos_0', 'mgw_outer_str_pos_1']
    xt_df = pd.DataFrame(xt)[topf]
    xtlist = xt_df.values.tolist()
    y_train = get_numeric_label(t.df).values
    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=5,random_state=0,min_samples_leaf=10,min_samples_split=10)
    m = rf.fit(xtlist, y_train)
    tree.export_graphviz(m.estimators_[5], out_file='tree.dot',
            feature_names = xt_df.columns,
            class_names = ['additive','cooperative'],
            rounded = True, proportion = False,
            precision = 2, filled = True)
    subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])

    # xt = xt_df.values.tolist()
    # rf.fit(xt,y_train)
    # pickle.dump(rf, open("dist_ori_flank_tt.sav", 'wb'))
