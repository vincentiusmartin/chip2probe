'''
Created on Oct 30, 2019

@author: vincentiusmartin
'''
import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL

//from chip2probe.modeler.oldtraining_gen.training import Training
from sklearn import ensemble, model_selection, metrics, tree
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
import pickle

import trainingdata.seqextractor as seqextractor
from trainingdata.dnashape import DNAShape

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
    plt.title('Average ROC curves using CART')
    leg = plt.legend(loc="lower right",title="AUC for each combination:")
    leg._legend_box.align = "left"
    plt.savefig(path)


def merge_listdict(ld1, ld2):
    """

    """
    if len(ld1) > 0 and len(ld2) > 0 and len(ld1) != len(ld2):
        print("Error: list lengths are not the same")
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

def make_linker_table(t, xlist):
    df = t.df[["sequence","label"]]
    to_appends = [df, pd.DataFrame(t.get_linker_list())]
    for xl in xlist:
        to_appends.append(pd.DataFrame(xl))
    df = pd.concat(to_appends, axis = 1)
    return df

# TODO FIX BASED ON THE CLASS
"""
def plot_average_all(train,shape,distances):
    for dist in distances:
        print("Plotting for dist %d" % dist)
        dist_path = "%s/d%s" % (shapepath,dist)
        # make a new data frame with only the distance on each iteration
        t2 = train.training.loc[train.training['distance'] == dist]
        train2 = trainingparser.TrainingParser(t2,motiflen=6)
        li = train2.get_labels_indexes()
        bsites = train2.get_bsites()
        shape = DNAShape(dist_path)

        for p in [0.05,0.1]: # 0.05,
            plot_path = "%s/shape-p=%.2f.pdf"%(dist_path,p)
            shape.plot_average(li,bsites,pthres=p,path=plot_path,plotlabel="Average DNA shape for d=%d,p=%.2f" % (dist,p))
"""

def plot_average_all(df,shapepath,distances,corelen):
    shape = DNAShape(shapepath)
    for dist in distances:
        print("Plotting for dist %d" % dist)
        dist_path = "%s/d%s" % (shapepath,dist)
        # make a new data frame with only the distance on each iteration
        df_dist = df.loc[dftrain['distance'] == dist]
        if df_dist.empty:
            continue
        labeled_idxs = df_dist.groupby("label").groups
        s1list,s2list = get_ordered_site_list(df_dist)
        seqdict = df_dist["sequence"].to_dict()
        for p in [0.05]: # 0.05,
            plot_path = "%s/shape-d%d-p=%.2f.pdf"%(shapepath,dist,p)
            shape.plot_average(labeled_idxs, s1list, s2list, seqdict, pthres=p,path=plot_path,plotlabel="Average DNA shape for d=%d,p=%.2f" % (dist,p))

def get_ordered_site_list(df):
    site1 = {}
    site2 = {}
    for idx,row in df.iterrows():
        if row["site_wk_pos"] > row["site_str_pos"]:
            site1[idx] = row["site_str_pos"]
            site2[idx] = row["site_wk_pos"]
        else:
            site1[idx] = row["site_wk_pos"]
            site2[idx] = row["site_str_pos"]
    return site1,site2

def stacked_bar_categories(df, x, y=["label"],plotname="stackedbar.png",avg=False,legend=True, ylabel=""):
    cat_df = df[[x]]
    #cat_df["orientation"] = cat_df["orientation"].map({'1': 'HT/TH', '2': 'HH', '3':'TT'})
    cat_df["label"] = df['label']

    #order = ['HT/TH', 'HH', 'TT']
    #.set_index('orientation').loc[order]

    group = [x] + y
    df2 = cat_df.groupby(group)['label'].count() # .unstack(x).fillna(0)
    if avg:
        df2 = df2.groupby(level=0).apply(lambda x: x / float(x.sum()))
    ax = df2.unstack(x).fillna(0).T.plot(kind='bar', stacked=True,legend=legend,rot=0)
    ax.set_ylabel(ylabel)
    plt.savefig(plotname)
    plt.clf()

def plot_auc(x_train_dict, df, plotname="auc.png"):
    #fpr_dict = {key:[] for key in x_train}
    tpr_dict = {key:[] for key in x_train_dict}
    auc_dict = {key:[] for key in x_train_dict}
    acc_dict = {key:[] for key in x_train_dict}

    fpr_lim=100
    base_fpr = np.linspace(0, 1, 101)[:fpr_lim+1]
    cv = model_selection.KFold(n_splits=10,shuffle=True)

    random_x = next(iter(x_train_dict.values()))[0]
    x_train = pd.DataFrame(random_x).values.tolist()
    y_train = get_numeric_label(df).values

    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=10,random_state=0)
    dt = tree.DecisionTreeClassifier(min_samples_split=27, min_samples_leaf=25, criterion="entropy")

    for train_idx,test_idx in cv.split(x_train):
        tprs = []
        aucs_val = []
        # compare with just distance
        for key in x_train_dict:
            # need to convert this with index, somehow cannot do
            # x_train[train_idx] for multi features
            xt = x_train_dict[key][0]
            xt = [[d[k] for k in d] for d in xt] # need to be a list of list
            model_name = x_train_dict[key][1]

            data_train = [xt[i] for i in train_idx]
            data_test = [xt[i] for i in test_idx]
            lbl_train = [y_train[i] for i in train_idx]
            lbl_test = [y_train[i] for i in test_idx]

            if model_name == "rf":
                model = rf.fit(data_train, lbl_train)
            else:
                model = dt.fit(data_train, lbl_train)
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

    # left append 0 in base fpr just so we start at 0 (we did the same for the tpr)s
    display_output(base_fpr, mean_tpr, mean_auc, path=plotname)

if __name__ == '__main__':
    #trainingpath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191004_coop-PBM_Ets1_v1_1st/training_data/training_overlap.tsv"
    trainingpath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/training_data/training_p01_adjusted.tsv"
    pd.set_option('display.max_columns', None)
    dforig = pd.read_csv(trainingpath, sep="\t")
    ori_one_hot = True
    feature_dist_type = "numerical"

    # only get cooperative and additive
    dforig = dforig[(dforig["label"] == "cooperative") | (dforig["label"] == "additive")]
    dftrain = dforig[~dforig['name'].str.contains(r'weak|dist')].reset_index(drop=True)# r'dist|weak
    tr1 = Training(dftrain, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])
    dftr = tr1.df

    x_ori = tr1.get_feature_orientation(["GGAA","GGAT"], one_hot = False)
    dftr["orientation"] = pd.DataFrame(x_ori)["ori"]
    df_ht = dftrain #[dftrain["orientation"] == "HT/TH"] #[dftrain["distance"] % 2 == 0] #[dftrain["orientation"] == "HT/TH"]
    dftr.to_csv("train1.tsv",sep="\t")

    t = Training(df_ht, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])
    #t.stacked_bar_categories("distance",avg=True)

    s_in = 5
    s_out = 4
    ds = DNAShape("/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/dnashape/training_p01_adjusted_reversed")

    # ========== GETTING FEATURES FROM THE DATA ==========

    x_linker_pos = t.get_linker_positional_feature(6, maxk=3)
    x_linker_shape = t.get_linker_positional_feature(6, dnashape=ds)

    x_dist_numeric = t.get_feature_distance(type="numerical")
    x_ori = t.get_feature_orientation(["GGAA","GGAT"], one_hot = ori_one_hot)
    x_gc = t.get_linker_GC_content()
    x_pref = t.get_feature_site_pref()
    smode = "strength" # or positional
    x_flank_in = t.get_feature_flank_core(k=3, seqin = s_in, site_mode=smode)
    x_flank_out = t.get_feature_flank_core(k=3, seqin = -s_out, site_mode=smode)
    x_shape_in = t.get_feature_flank_shapes(ds, seqin = s_in, site_mode=smode)
    x_shape_out = t.get_feature_flank_shapes(ds, seqin = -s_out, site_mode=smode)

    x_mid_shape_mean = t.get_middle_avgshape_feature([1,3,5], ds ,maxk=2, action="mean")
    x_mid_shape_max = t.get_middle_avgshape_feature([1,3,5], ds ,maxk=2, action="max")
    x_mid_shape_min = t.get_middle_avgshape_feature([1,3,5], ds ,maxk=2, action="min")

    xtr = []
    for x in [x_dist_numeric, x_ori, x_flank_in, x_flank_out]: #
        xtr = merge_listdict(xtr, x)
    x_df = pd.DataFrame(xtr)
    x_df.to_csv("features.csv", index=False)

    # ========== CREATING THE RF OBJECT  ==========
    x_train = pd.DataFrame(xtr).values.tolist()
    y_train = get_numeric_label(t.df).values

    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=10,random_state=0)
    dt = tree.DecisionTreeClassifier(min_samples_split=27, min_samples_leaf=25, criterion="entropy")
    #y_pred = model.predict(x_train)

    # ========== GET TOP N FEATURES  ==========
    model = rf.fit(x_train, y_train)
    feature_importances = pd.DataFrame(rf.feature_importances_,
                                   index = x_df.columns,
                                   columns=['importance']).sort_values('importance', ascending=False)
    imp = list(feature_importances.index[:10])
    print("Top 10 feature importance list " + str(imp))
    x_df_imp = x_df[imp].to_dict('records') # we only use the most important features
    xtr_imp =  merge_listdict([], x_df_imp)

    #x_train_dict = {"top5": x_df_imp.values.tolist()} #, "dist-numeric":x_df[["dist-numeric"]].values.tolist()}

    x1 = merge_listdict([], x_dist_numeric)

    xt2 = []
    for x in [x_dist_numeric, x_linker_shape]:
        xt2 = merge_listdict(xt2, x)

    xt3 = []
    for x in [x_dist_numeric, x_linker_pos]:
        xt3 = merge_listdict(xt3, x)
    """
    for x in [x_shape_out, x_flank_in, x_shape_in, x_flank_out]: #[x_dist_numeric,x_ori,x_link1,x_link2,x_link3,x_gc,x_pref]: #  x_pref,
        x2 = merge_listdict(x2, x)

    x3 = []
    for x in [x_dist_numeric, x_ori, x_pref]: #[x_dist_numeric,x_ori,x_link1,x_link2,x_link3,x_gc,x_pref]: #  x_pref,
        x3 = merge_listdict(x3, x)
    """

    x_imp = merge_listdict([], x_df_imp)
    #x_train_dict  = {"distance,strength,orientation":[x6,"dt"], "distance,orientation":[x5,"dt"], "distance,strength":[x4,"dt"], "site strength":[x3,"dt"], "orientation":[x2,"dt"], "distance":[x1,"dt"]}
    x_train_dict = {"distance":[x1,"dt"], "all":[xtr,"dt"], "top10": [xtr_imp,"rf"]}

    # make_model
    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=10,random_state=0)
    dt = tree.DecisionTreeClassifier(min_samples_split=27, min_samples_leaf=25, criterion="entropy")
    
    xt = [[d[k] for k in d] for d in xtr]
    yt = get_numeric_label(t.df).values
    dt.fit(xt,yt)
    pickle.dump(rf, open("model1_all_all_dt.sav", 'wb'))
    plot_auc(x_train_dict,t.df)

    # ========== TREE DRAWING FROM A MODEL  ==========
    """
    x_dict_tree = xtr_imp
    x_df_tree = pd.DataFrame(x_dict_tree)
    x_list_tree = [[d[k] for k in x_df_tree.columns] for d in x_dict_tree]
    model = dt.fit(x_list_tree, y_train) # make new model from the top features
    # take a tree, let's say tree #5
    estimator = model
    tree.export_graphviz(estimator, out_file='tree.dot',
            feature_names = x_df_tree.columns,
            class_names = ['additive','cooperative'],
            rounded = True, proportion = False,
            precision = 2, filled = True)
    subprocess.call(['dot', '-Tpdf', 'tree.dot', '-o', 'tree.pdf', '-Gdpi=600'])

    #print("Accuracy on test: %.f" % accuracy_score(y_true, y_pred))


    # ========== Using this result to train on different dataset  ==========

    testpath = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/training_data/training_with_coop_anti.tsv"
    dftest = pd.read_csv(testpath, sep="\t")
    dftest = dftest[(dftest["label"] == "cooperative") | (dftest["label"] == "additive")]

    test = Training(dftest, corelen=4)
    xtest_dist = test.get_feature_distance(type=feature_dist_type)
    xtest_link1 = test.get_feature_linker_composition(1)
    xtest_link2 = test.get_feature_linker_composition(2)
    xtest_link3 = test.get_feature_linker_composition(3)
    #xtest_gc = test.get_linker_GC_content()
    xtest_ori = test.get_feature_orientation(["GGAA","GGAT"], one_hot = ori_one_hot)
    xtest_pref = test.get_feature_site_pref()

    x_test = []
    for x in [xtest_dist,xtest_link1,xtest_link2,xtest_link3,xtest_ori,xtest_pref,xtest_gc]:
        x_test = merge_listdict(x_test, x)
    x_test = pd.DataFrame(x_test)[imp].values.tolist()
    y_true = get_numeric_label(t.df).values

    y_pred = model.predict(x_test)
    lpred = ["cooperative" if p == 1 else "additive" for p in y_pred]
    dfpred = dftest[["name"]]
    dfpred["pred"] = lpred
    dfpred.to_csv("aa.csv")
    """
