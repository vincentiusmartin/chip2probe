'''
Created on Oct 30, 2019

@author: vincentiusmartin
@editedby: Farica Zhuang
@description: all plot modules for a cooperative training data
'''


import matplotlib.pyplot as plt
import seaborn as sns
import math
import itertools
import numpy as np
import scipy
from decimal import Decimal
import os
from sklearn import model_selection, metrics
import chip2probe.util.stats as st

def plot_stacked_categories(df, x, y="label", path="stackedbar.png",
                           avg=False, legend=True):
    """
    Plot a stacked bar graph for each label.

    The groups are in the x-axis and counts are in the y-axis.
    Each bar is separated into different colors based on the labels.

    Args:
        x : Feature to group.
        y : The y axis, default is label
        path : Where to save the plot
        avg: If True then calculate the ratio between value in x Axis.
        legend: If True then show legend

    Returns:
        (NA)

    Example:
    """

    cat_df = df[[x] + [y]]
    group = [x] + [y]
    df2 = cat_df.groupby(group)[y].count() # .unstack(x).fillna(0)
    if avg:
        df2 = df2.groupby(level=0).apply(lambda x: x / float(x.sum()))
    ax = df2.unstack(x).fillna(0).T.plot(kind='bar', stacked=True,
                                         legend=legend, rot=0)
    ylabel = "ratio" if avg else "count"
    ax.set_ylabel(ylabel)
    plt.savefig(path)
    plt.clf()

def plot_box_categories(df, by=["label"], incols="default", path="boxplot.png", alternative="greater"):
    """
    Make boxplot.

    desc

    Args:


    Returns:
        (NA)

    Example:
        plot_box_categories(df, incols=["distance", "site_str_score", "site_wk_score"])
    """
    if incols == "default":
        cols = list(set(df.columns) - set(by))
    else:
        cols = list(incols)
    sns.set_style("whitegrid")
    numcol = 4
    numrow = math.ceil(len(cols) / numcol)
    # to make axis with different y-scale
    fig, ax = plt.subplots(numrow, numcol, figsize=(14, 5))
    plt.subplots_adjust(hspace =0.4, wspace=0.6)
    grouped = df.groupby(by=by)
    # need to sort to keep the order consistent
    cols.sort()
    logstr = ""
    for i in range(len(cols)):
        colname = cols[i]
        cur_group = {elm[0]:list(elm[1]) for elm in grouped[colname]}
        labels, data = [*zip(*cur_group.items())]
        cur_ax = ax.flatten()[i]
        sns.boxplot(data=data, width=.5, ax=cur_ax)
        # sns.stripplot(data=data, jitter=True, ax=cur_ax, size=1, color='k')
        cur_ax.set_xticklabels(labels)
        cur_ax.set_title(colname)
        #boxplot = self.df.boxplot(column=colname,by=by,ax=ax.flatten()[i],return_type="dict") # pandas version

        # only for group of size 1 for now:
        combs = list(itertools.combinations(range(len(data)), 2))
        ylim = float("-inf")
        for comb in combs:
            x1, x2 = comb[0],comb[1]
            max1 = max(data[x1])
            max2 = max(data[x2])
            mval = max1 if max1 > max2 else max2
            p_gr = st.wilcox(data[x1],data[x2],alternative=alternative)
            hfrac = -1 if mval < 0 else 1
            hfactor = (max2 - max1)**1.5
            pline_h = mval * 0.05
            pline_pos = mval * 0.2
            y, h, col = df[colname].max() + hfactor * pline_pos, pline_h, 'k'
            logstr += "%s, %s > %s: %4E\n" % (colname, labels[0], labels[1], Decimal(p_gr))
            pstr = "%.2E" % Decimal(p_gr) if p_gr < 0.001 else "%.4f" % p_gr
            cur_ax.plot([x1, x1, x2, x2], [1.01*y, y+h, y+h, 1.01*y], lw=1, c=col)
            cur_ax.text((x1+x2)*.5, y + h, "p = %s"%(pstr), ha='center', va='bottom', color="red")

            adj_ylim = (y+h) * 1.1
            if adj_ylim > ylim: ylim = adj_ylim
        cur_ax.set_ylim(top=ylim)
    for d in range(len(cols),numrow*numcol):
        fig.delaxes(ax.flatten()[d])
    plt.savefig(path)
    plt.clf()
    # also write the log file
    basepath = os.path.splitext(path)[0]
    with open("%s.log" % basepath, 'w') as f:
        f.write(logstr)

def display_output(xy, score_dict, path, title, score_type="auc", varyline=False):
    """
    Display output

    This plots the average ROC curve of all the classifiers in a single plot

    Args:
        xy: a dictionary of dictionary with the plot title as key. The subdictionary
           contains 'x' and 'y' where each is a list of values in the plot.
        score_dict: a dictionary of auc/pr value, the keys in xy should be the same
           with keys in score_dict
        path:
        title:
        score_type:
        varyline: if True generates different type of line depending # features
    Returns:
        (NA)
    """
    plt.clf() # first, clear the canvas
    plt.plot([0, 1], [0, 1], color="gray", alpha=0.5, lw=0.3)#linestyle="--",
    for key in  xy:
        score = score_dict[key]
        if varyline:
            ln = key.split(",")
            if len(ln) >= 3:
                lw, ls = 2, "-"
            elif len(ln) == 2:
                lw, ls = 1, "-"
            else:
                lw, ls = 1, "--"
        else:
            lw, ls = 1, "-"

        plt.plot(xy[key]['x'], xy[key]['y'], lw=lw, linestyle=ls , label='%s: %.2f' % (key,score))

        # Show the ROC curves for all classifiers on the same plot
        if score_type == "pr":
            plt.xlabel('Recall')
            plt.ylabel('Precision')
        else:
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')

    plt.title(title)
    leg = plt.legend(loc="lower right",title="%s for each combination:"%str.upper(score_type))
    leg._legend_box.align = "left"
    plt.savefig(path, dpi=600)


def plot_model_metrics(modeldict, plotname="auc.png",
                 title="Average ROC Curves", cvfold=10, score_type="auc",
                 interp=False, writelog=True, varyline = False):
    """
    Make boxplot.

    x== y==

    Args:
        -modeldict: a dictionary with the model name as the key and  features,
                model as value. This can be obtained using BestModel class
        -plotname: where to save the plot
        -title: plot title
        -cvfold: cross validation fold
        -score_type: auc/pr
        -interp: whether to interpolate (smooth the curve)
        -writelog: write metrics log to file
        -varyline: if True generates different type of line depending # features
    Returns:
        (NA)
    """
    if interp and score_type == "pr":
        raise ValueError("Could not interpolate with precision recall")

    yprob_dict = {key:[] for key in modeldict}
    ytrue_dict = {key:[] for key in modeldict}
    ypred_dict = {key:[] for key in modeldict}

    cv = model_selection.KFold(n_splits=cvfold,shuffle=True)

    random_x = next(iter(modeldict.values()))[0]
    xt_all = random_x.values.tolist()
    labels = random_x['label'].values

    for train_idx,test_idx in cv.split(xt_all):
        aucs_val = []
        for key in modeldict:
            # need to convert this with index, somehow couldn't do
            # x_train[train_idx] for multi features
            xt_df = modeldict[key][0]
            xt = np.array(xt_df.loc[:,xt_df.columns != 'label'].values.tolist()) # need to be a list of list

            x_train, x_test = xt[train_idx], xt[test_idx]
            y_train, y_test = labels[train_idx], labels[test_idx]

            model = modeldict[key][1]
            model = model.fit(x_train, y_train)
            # if key == "Lasso" or key == "Elastic Net":
            #     y_score = model.predict(x_test)
            #     y_label = [1 if i >= 0.5 else 0 for i in y_score]
            # else:
            y_score = model.predict_proba(x_test)[:,1]
            y_label = model.predict(x_test)

            ytrue_dict[key].append(y_test)
            yprob_dict[key].append(y_score)
            ypred_dict[key].append(y_label)

    if interp: #interpolate the score
        base_x = np.linspace(0, 1, 101)
        scoreval = {}
        mets = {}
        for k in ytrue_dict:
            scores = []
            ylist = []
            for i in range(0,len(ytrue_dict[k])):
                # we calculate the auc and pr score before interpolation for better representation
                if score_type == "auc":
                    x, y, _ = metrics.roc_curve(ytrue_dict[k][i], yprob_dict[k][i]) #fpr, tpr, threshold
                    scores.append(metrics.auc(x, y))
                elif score_type == "pr":
                    y, x, _ = metrics.precision_recall_curve(ytrue_dict[k][i], yprob_dict[k][i]) # precision, recall, threshold
                    scores.append(metrics.average_precision_score(x, y))
                y = scipy.interp(base_x, x, y) # add interp points to the plotting
                ylist.append(y) # insert 0 so we start the plot from the bottom left
            avgscr = np.array(ylist).mean(axis=0)
            mets[k] = {
                       'x': np.insert(base_x,0,0),
                       'y': np.insert(avgscr,0,0),
                       }
            scoreval[k] = np.mean(scores)
    else:
        ytrue_dict = {k:np.concatenate(ytrue_dict[k]) for k in ytrue_dict}
        yprob_dict = {k:np.concatenate(yprob_dict[k]) for k in yprob_dict}
        if score_type == "auc":
            mets = {k: metrics.roc_curve(ytrue_dict[k], yprob_dict[k]) for k in ytrue_dict} #fpr, tpr, threshold
            for k in mets:
                mets[k] = {"x":mets[k][0], "y":mets[k][1]} #"fpr, tpr"
            scoreval = {k:metrics.auc(mets[k]['x'],mets[k]['y']) for k in ytrue_dict} #{k:np.array(auc_dict[k]).mean(axis=0) for k in auc_dict}
        elif score_type == "pr":
            mets = {k: metrics.precision_recall_curve(ytrue_dict[k], yprob_dict[k]) for k in ytrue_dict} # precision, recall, threshold
            for k in mets:
                mets[k] = {"x":mets[k][1], "y":mets[k][0]} #recall, precision
            scoreval = {k:metrics.average_precision_score(ytrue_dict[k],yprob_dict[k]) for k in ytrue_dict}

    ypred_dict = {k:np.concatenate(ypred_dict[k]) for k in ypred_dict}
    if interp: # because we don't have this during interp:
        ytrue_dict = {k:np.concatenate(ytrue_dict[k]) for k in ytrue_dict}
    acc = {k:metrics.accuracy_score(ytrue_dict[k],ypred_dict[k]) for k in ytrue_dict}
    confmat = {k:metrics.confusion_matrix(ytrue_dict[k],ypred_dict[k]).ravel() for k in ytrue_dict}

    # TODO: make better represetation of this
    logstr = (
            "Mean accuracy: %s\n" +
            "Mean %s: %s\n" +
            "Confusion_matrix (tn,fp,fn,tp): %s"
            ) % (str(acc), score_type, str(scoreval), str(confmat))
    if writelog:
        logpath = "%s.log" % os.path.splitext(plotname)[0]
        with open(logpath, 'w') as f:
            f.write(logstr)
    display_output(mets, scoreval, path=plotname, title=title, score_type=score_type, varyline=varyline)
