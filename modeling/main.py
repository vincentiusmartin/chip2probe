'''
Created on Oct 30, 2019

@author: vincentiusmartin
'''
import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL

from trainingdata.training import Training
from sklearn import ensemble, model_selection, metrics, tree
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score

import subprocess

def get_numeric_label(training):
    # hard coded but change add to anti coop / additive when needed
    train = training['label'].map({'cooperative': 1, 'additive': 0})
    return train

def display_output(fpr_list, tpr_dict, auc_dict, path):
    """
        This plots the average ROC curve of all the classifiers in a single plot
    """
    plt.clf() # first, clear the canvas
    for key in tpr_dict:
        tpr_list = tpr_dict[key]
        auc = auc_dict[key]

        plt.plot([0, 1], [0, 1], linestyle="--", color="red", alpha=0.1)
        #for key in fpr_dict:
        plt.plot(fpr_list, tpr_list, lw=2, alpha=0.4, label='%s, AUC %f' % (key,auc))

        # Show the ROC curves for all classifiers on the same plot
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
    plt.title('Average ROC Curves for All Classifiers')
    plt.legend(loc="lower right")
    plt.savefig(path)


def merge_listdict(ld1, ld2):
    if len(ld1) > 0 and len(ld2) > 0 and len(ld1) != len(ld2):
        print("Error:list length is not the same")
        return 0
    if len(ld1) == 0:
        return ld2
    elif len(ld2) == 0:
        return ld1
    else:
        ld_ret = []
        for i in range(0,len(ld1)):
            ld_ret.append({**ld1[i], **ld2[i]})
        return ld_ret

def get_custom_df(df,keystr):
    # first get df with dist custom
    df_cust = df[df['name'].str.contains(keystr)]

    all_names = set(map(lambda x: "_".join(x.split("_")[:-2]),df_cust["name"]))
    df_wt = df[df['name'].isin(all_names)]
    if keystr == "dist":
        df_wt["name"] = df_wt.apply(lambda row: row["name"] + "_dist_" + str(row["distance"]), axis=1)
    else:
        df_wt["name"] = df_wt.apply(lambda row: "%s_%s" % (row["name"],keystr), axis=1)
    df_cust_all = pd.concat([df_cust, df_wt], axis=0)
    return df_cust_all

if __name__ == '__main__':
    #trainingpath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191004_coop-PBM_Ets1_v1_1st/training_data/training_overlap.tsv"
    trainingpath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/training_data/training_overlap.tsv"
    pd.set_option('display.max_columns', None)
    dforig = pd.read_csv(trainingpath, sep="\t")
    dforig = dforig[(dforig["label"] == "cooperative") | (dforig["label"] == "additive")]

    # only get cooperative and additive
    dftrain = dforig[~dforig['name'].str.contains(r'dist|weak')]
    #
    #dftrain = get_custom_df(dforig,"dist|weak")

    t = Training(dftrain, corelen=4)
    #t.plot_summary()

    # ========== GETTING FEATURES FROM THE DATA ==========
    x_dist = t.get_feature_distance(type="numerical")
    x_link1 = t.get_feature_linker_composition(1)
    x_link2 = t.get_feature_linker_composition(2)
    x_link3 = t.get_feature_linker_composition(3)
    x_ori = t.get_feature_orientation(["GGAA","GGAT"])
    print(x_ori)
    x_pref = t.get_feature_site_pref()

    x_train = []
    for x in [x_dist,x_link1,x_link2,x_link3,x_ori,x_pref]: # , x_link1, x_link2, x_link3, x_orix`
        x_train = merge_listdict(x_train, x)
    x_df = pd.DataFrame(x_train)

    # ========== CREATING THE RF OBJECT  ==========
    x_train = pd.DataFrame(x_train).values.tolist()
    y_train = get_numeric_label(t.df).values

    rf = ensemble.RandomForestClassifier(n_estimators=100, max_depth=10,random_state=0)
    #y_pred = model.predict(x_train)

    # ========== GET TOP N FEATURES  ==========
    model = rf.fit(x_train, y_train)
    feature_importances = pd.DataFrame(rf.feature_importances_,
                                   index = x_df.columns,
                                   columns=['importance']).sort_values('importance', ascending=False)
    imp = list(feature_importances.index[:5])
    print("Top 10 feature importance list " + str(imp))
    x_df_imp = x_df[imp] # we only use the most important features
    x_train_dict = {"top5": x_df_imp.values.tolist(), "dist-numeric":x_df[["dist-numeric"]].values.tolist()}

    # ========== TREE DRAWING FROM A MODEL  ==========
    model = rf.fit(x_train_dict["top5"], y_train) # make new model from the top features

    # take a tree, let's say tree #5
    estimator = rf.estimators_[5]
    tree.export_graphviz(estimator, out_file='tree.dot',
            feature_names = x_df_imp.columns,
            class_names = ['cooperative', 'additive'],
            rounded = True, proportion = False,
            precision = 2, filled = True)
    subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])
    # ========== Using this result to train on different dataset  ==========
    testpath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/training_data/training_with_coop_anti.tsv"
    dftest = pd.read_csv(testpath, sep="\t")
    dftest = dftest[(dftest["label"] == "cooperative") | (dftest["label"] == "additive")]

    test = Training(dftest, corelen=4)
    xtest_dist = test.get_feature_distance(type="numerical")
    xtest_link1 = test.get_feature_linker_composition(1)
    xtest_link2 = test.get_feature_linker_composition(2)
    xtest_link3 = test.get_feature_linker_composition(3)
    xtest_ori = test.get_feature_orientation(["GGAA","GGAT"])
    xtest_pref = test.get_feature_site_pref()

    x_test = []
    for x in [xtest_dist,xtest_link1,xtest_link2,xtest_link3,xtest_ori,xtest_pref]:
        x_test = merge_listdict(x_test, x)
    x_test = pd.DataFrame(x_test)[imp].values.tolist()
    y_true = get_numeric_label(t.df).values

    y_pred = model.predict(x_test)
    lpred = ["cooperative" if p == 1 else "additive" for p in y_pred]
    dfpred = dftest[["name"]]
    print(len(lpred))
    dfpred["pred"] = lpred
    dfpred.to_csv("aa.csv")
    #print("Accuracy on test: %.f" % accuracy_score(y_true, y_pred))


    # ========== MAKING AUC USING THE TOP FEATURES  ==========

    #fpr_dict = {key:[] for key in x_train}
    tpr_dict = {key:[] for key in x_train_dict}
    auc_dict = {key:[] for key in x_train_dict}

    fpr_lim=100
    base_fpr = np.linspace(0, 1, 101)[:fpr_lim+1]
    cv = model_selection.KFold(n_splits=10,shuffle=True)

    for train_idx,test_idx in cv.split(x_train,y_train):
        tprs = []
        aucs_val = []
        # compare with just distance
        for key in x_train_dict:
            # need to convert this with index, somehow cannot do
            # x_train[train_idx] for multi features
            xt = x_train_dict[key]
            data_train = [xt[i] for i in train_idx]
            data_test = [xt[i] for i in test_idx]
            lbl_train = [y_train[i] for i in train_idx]
            lbl_test = [y_train[i] for i in test_idx]

            model = rf.fit(data_train, lbl_train)
            y_score = model.predict_proba(data_test)
            fpr, tpr, _ = metrics.roc_curve(lbl_test, y_score[:, 1])
            #auc = metrics.roc_auc_score(lbl_test, y_score[:,1])
            tpr = scipy.interp(base_fpr, fpr, tpr) # add points to the plotting
            res_auc = metrics.auc(base_fpr, tpr)
            tpr_dict[key].append(np.insert(tpr,0,0)) # insert 0 so we start the plot from the bottom left
            auc_dict[key].append(res_auc)

    mean_tpr = {k:np.array(tpr_dict[k]).mean(axis=0) for k in tpr_dict}
    mean_auc = {k:np.array(auc_dict[k]).mean(axis=0) for k in auc_dict}

    # left append 0 in base fpr just so we start at 0 (we did the same for the tpr)
    base_fpr = np.insert(base_fpr,0,0)
    display_output(base_fpr, mean_tpr, mean_auc, path="here.png")
