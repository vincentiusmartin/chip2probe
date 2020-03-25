import sys
#sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL
sys.path.append("/Users/faricazjj/Desktop/homotf/chip2probe")
from trainingdata.training import Training
import util.util as util
from best_model import BestModel

from sklearn import ensemble, model_selection, metrics, tree
import pandas as pd
import numpy as np
import scipy
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt

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

def plot_auc(x_train_dict, df, plotname="auc.png"):
    tpr_dict = {key:[] for key in x_train_dict}
    auc_dict = {key:[] for key in x_train_dict}
    acc_dict = {key:[] for key in x_train_dict}

    fpr_lim=100
    base_fpr = np.linspace(0, 1, 101)[:fpr_lim+1]
    cv = model_selection.KFold(n_splits=10,shuffle=True)

    random_x = next(iter(x_train_dict.values()))[0]
    x_train = pd.DataFrame(random_x).values.tolist()

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

            model = x_train_dict[key][1]

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

    # left append 0 in base fpr just so we start at 0 (we did the same for the tpr)
    base_fpr = np.insert(base_fpr,0,0)
    display_output(base_fpr, mean_tpr, mean_auc, path=plotname)

def get_top_n(n, xdict, ytrain, rf):
    x_df = pd.DataFrame(xdict)
    xtrain = x_df.values.tolist()
    model = rf.fit(xtrain, ytrain)
    feature_importances = pd.DataFrame(rf.feature_importances_,
                                   index = x_df.columns,
                                   columns=['importance']).sort_values('importance', ascending=False)
    imp = list(feature_importances.index[:n])
    print("Top 10 feature importance list " + str(imp))
    x_df_imp = x_df[imp].to_dict('records') # we only use the most important features
    return  util.merge_listdict([], x_df_imp)

if __name__ == "__main__":
    #trainingpath = "train1.tsv"
    trainingpath = "trainingdata/training_new.csv"

    df = pd.read_csv(trainingpath, sep=",")

    t = Training(df, corelen=4).flip_one_face_orientation(["GGAA","GGAT"])
    y_train = get_numeric_label(t.df).values

    rf = ensemble.RandomForestClassifier(n_estimators=500, max_depth=10,random_state=0)
    dt = tree.DecisionTreeClassifier(min_samples_split=27, min_samples_leaf=25, criterion="entropy")
    xtr = {"distance": [
                    t.get_feature_all({
                        "distance":{"type":"numerical"}
                    }),
                    dt
                ],
            "all": [
                    t.get_feature_all({
                        "distance":{"type":"numerical"},
                        "flankseq": {"k":3, "seqin":4, "smode":"strength"},
                        "flankseq": {"k":3, "seqin":-3, "smode":"strength"}
                    }),
                    rf
                ]
        }
    xtr["topn"] = [get_top_n(10, xtr["all"][0], y_train, rf), rf]

    plot_auc(xtr, y_train, "auc.png")

    """
    param_dict = {
    				'n_estimators': [i for i in range(2,21)],
    				'max_depth': [i for i in range(100,2001,100)]
    			}
    rf = BestModel(clf="RF", param_dict=param_dict, topn=10, train_data=train_data).run_all()
    """