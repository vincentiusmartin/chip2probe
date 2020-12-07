import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import statsmodels.stats.multitest as sm

import chip2probe.util.stats_r as st

def assign_fdrcor_class(p, prevlbl, pcut=0.05):
    if prevlbl == "fail_cutoff":
        return 'fail_cutoff'
    elif prevlbl == 'anticoop' and p < pcut:
        return 'anticoop'
    elif prevlbl == 'cooperative' and p < pcut:
        return 'cooperative'
    else:
        return 'additive'

def create_cooplbl(indivsum, twosites, pcutoff = 0.05):
    """
    """
    p_coop = st.wilcox(twosites, indivsum, "greater")
    p_anti = st.wilcox(twosites, indivsum, "less")
    if p_coop < pcutoff:
        return "cooperative", p_coop
    elif p_anti < pcutoff:
        return "anticoop", p_anti
    else:
        return 'additive', p_coop

def label_replicas_permutation(indiv, two, arrdf, cutoff=0, oricol="ori", namecol="Name", affcol="affinity", typecol="type", pcut = 0.05):
    median_dict = arrdf.groupby([namecol, oricol, typecol])[affcol].median().to_dict()
    labeled_dict = {}
    for ori in list(indiv.keys()):
        orilbls = []
        for k in indiv[ori]:
            rowdict = {}
            if median_dict[(k,ori,'wt')] < cutoff or median_dict[(k,ori,'m3')] > cutoff:
                rowdict['label'], rowdict['p'] = "fail_cutoff", 1
            else:
                rowdict['label'], rowdict['p'] =  create_cooplbl(indiv[ori][k], two[ori][k], pcut)
            rowdict['indiv_median'] = np.median(indiv[ori][k])
            rowdict['two_median'] = np.median(two[ori][k])
            rowdict['Name'] = k
            orilbls.append(rowdict)
        labeled_dict[ori] = pd.DataFrame(orilbls)
        labeled_dict[ori].to_csv("%s.csv"%ori,index=False)
        labeled_dict[ori]['p'] = sm.fdrcorrection(labeled_dict[ori]['p'])[1]
        labeled_dict[ori]['label'] = labeled_dict[ori].apply(lambda row: assign_fdrcor_class(row['p'],row['label'],pcut),axis=1)
    return labeled_dict

def make_replicas_permutation(df, oricol="ori", namecol="Name", affcol="affinity", typecol="type", repcol="rep"):
    """
    Make permutation across replicas

    Args:
        type: wt,m1,m2,m3
    """
    oris = df[oricol].unique()
    indivsum = {ori:{} for ori in oris}
    twosites = {ori:{} for ori in oris}
    for ori in oris:
        df_ori = df[df["ori"] == ori]
        probes_dict = dict(tuple(df_ori.groupby(namecol)))
        for key, curdf in probes_dict.items():
            # maybe do this just for the first time?
            if set(curdf[typecol].unique()) != {"wt", "m1", "m2", "m3"}:# skip if incomplete
                continue
            g = curdf.set_index([repcol,typecol]).to_dict()[affcol]
            reps = curdf.groupby(typecol)[repcol].apply(list).to_dict()

            indivreps = list(itertools.product(*[reps[t] for t in ['m1','m2','m3']]))
            indivsum_list = []
            for r in indivreps:
                i_s = g[(r[0], 'm1')] + g[(r[1], 'm2')] - g[(r[2], 'm3')]
                indivsum_list.append(i_s)
            indivsum[ori][key] = indivsum_list

            tworeps = list(itertools.product(*[reps[t] for t in ['wt']]))
            twosites_list = []
            for r in tworeps:
                t_s = g[(r[0], 'wt')]
                twosites_list.append(t_s)
            twosites[ori][key] = twosites_list
    return indivsum, twosites

# --------

def plot_chamber_corr(dfx, dfy, xlab="Chamber 1", ylab="Chamber 2",
                namecol="Name", valcol="Alexa488Adjusted", extrajoincols=[],
                title="", cutoff="",
                median=False, log=False,
                shownames=False, path=""):
    """
    Get correlation between 2 chambers

    Args:
        dfx: data frame from chamber 1
        dfy: data frame from chamber 2
        namecol: column name for the probe names
        valcol: column name for the array intensity
        extrajoincols: column name for the extras such as orientation
        title: plot title
        median: get median between replicates
        log: return result in log scale
        shownames: show names of the probe on the plot, useful for debugging
        path: if empty then show the plot
    Return:
        R^2
    """
    # correlation plot of negative controls in chamber 1 vs chamber 2

    joincols = [namecol] + extrajoincols
    dfcombined = dfx.merge(dfy, on=joincols)[joincols + ["%s_x"%valcol, "%s_y"%valcol]]

    if median:
        dfcombined = dfcombined.groupby([namecol] + extrajoincols)[["%s_x"%valcol, "%s_y"%valcol]].median().reset_index()
    if log:
        dfcombined["%s_x"%valcol] = np.log(dfcombined["%s_x"%valcol])
        dfcombined["%s_y"%valcol] = np.log(dfcombined["%s_y"%valcol])

    x = dfcombined["%s_x"%valcol].values
    y = dfcombined["%s_y"%valcol].values
    # x = [i/10000 for i in x]
    # y= [i/10000 for i in y]

    f = plt.figure()
    ax = f.add_subplot(111)
    # plot diagonal
    plt.plot([min(min(x),min(y)),max(max(x),max(y))], [min(min(x),min(y)),max(max(x),max(y))], color='blue', label='diagonal') # , linewidth=0.5

    if cutoff:
        c = np.log(cutoff) if log else cutoff
        plt.axhline(cutoff, color='black', linewidth=0.5, linestyle=":")
        plt.axvline(cutoff, color='black', linewidth=0.5, linestyle=":")

    slope, intercept = np.polyfit(x, y, 1)
    # Create a list of values in the best fit line
    abline_values = [slope * i + intercept for i in x]

    # plot best fit line
    plt.plot(x, abline_values, color='black', label='best fit')
    r_squared = np.corrcoef(x,y)[0,1]**2
    plt.scatter(x, y, s=1) # color="#FFA07A"

    if shownames:
        names = dfcombined[namecol].values
        for i in range(len(names)):
            plt.text(x[i], y[i], names[i], fontsize=3)

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))
    plt.title(title)
    plt.text(0.02, 0.94, 'R² = %0.2f' % r_squared, color='red', fontsize=12, transform = ax.transAxes)
    sign = '+' if intercept > 0 else '-'
    plt.text(0.02, 0.9, 'best fit: y=%0.2fx %s %0.2f' % (slope,sign,abs(intercept)), color='red', transform = ax.transAxes, fontsize=12)
    plt.legend(loc='lower right', prop={'size': 11})

    if path:
        plt.savefig(path)
    else:
        plt.show()
    plt.clf()
    return r_squared

def plot_classified_labels(df, path="", col1="Alexa488Adjusted_x", col2="Alexa488Adjusted_y",
                           xlab="Chamber1", ylab="Chamber2", log=True, title="", axes=None,
                           shownames=False, namecol="Name", labelcol="label", plotnonsignif=True,
                           labelnames=["cooperative","additive","anticoop"]):
    """
    Desc

    Args:
        df: with "Name", "Intensity_one", "Intensity_two", "label".
            Label -> cooperative, additive, anticoop, fail_cutoff
            shownames: show names of the probe on the plot, useful for debugging. Names are obtained from namecol.

    Return:

    """
    permitted_labels = list(labelnames)
    if plotnonsignif:
        permitted_labels.append("fail_cutoff")
    newdf = df[df[labelcol].isin(permitted_labels)]

    ax = axes if axes else plt.axes()

    # red/firebrick, mistyrose, blue, skyblue
    lblclr = [("fail_cutoff","yellow", 0.7), (labelnames[1],"gray",1), (labelnames[0],"blue",1), (labelnames[2],"red",1)]
    if log:
        newdf[col1] = np.log(df[col1])
        newdf[col2] = np.log(df[col2])

    for lc in lblclr:
        if not newdf[newdf[labelcol]==lc[0]].empty:
            x = newdf[newdf[labelcol]==lc[0]][col1].values
            y = newdf[newdf[labelcol]==lc[0]][col2].values
            label = lc[0] if lc[0] != "fail_cutoff" else None
            ax.scatter(x, y, color=lc[1], s=3, alpha=lc[2], label=label)
            if shownames:
                names = newdf[newdf['label']==lc[0]][namecol].values
                for i in range(len(names)):
                    plt.text(x[i], y[i], names[i], fontsize=5)

    ax.set_xlabel(xlab, fontsize=12)
    ax.set_ylabel(ylab, fontsize=12)
    x,y = newdf[col1].values, newdf[col2].values
    ax.plot([min(min(x),min(y)),max(max(x),max(y))], [min(min(x),min(y)),max(max(x),max(y))], color='black')
    ax.legend(loc="lower right", markerscale=5, fontsize=10)
    ax.set_title(title, fontsize=15)

    if not axes:
        if not path:
            plt.show()
        else:
            plt.savefig(path)
        plt.clf()

def scatter_boxplot(input_df, cols=[], log=False,
                    ax=None, path = "box.png", repcol="rep",
                    title=""):
    """
    Scatter boxplot of intensity among replicates

    Args:
        input_df: df input
        path: where to save the plot
        cols: columns with intensity
    Return:

    """
    df = input_df[cols] if cols else input_df.copy()

    if not ax:
        df.boxplot(grid=False)
        p = plt
    else:
        p = df.boxplot(grid=False, ax=ax)

    for i in range(len(cols)):
        y = df[cols[i]]
        x = np.random.normal(i+1, 0.04, size=len(y))
        p.plot(x, y, 'r.', alpha=0.2)

    if ax:
        p.set_title(title)
    else:
        p.title(title)
        p.savefig(path)
        p.clf() # clear canvas

def scatter_boxplot_col(input_df, coltype, colval, plotorder=[], colororder=[], ax=None, path = "box.png", title=""):
    """
    """
    plotlist = []
    for po in plotorder:
        plotlist.append(input_df[input_df[coltype] == po][colval].tolist())

    idxs = list(range(1,len(plotorder)+1))
    if ax:
        ax.boxplot(plotlist)
        ax.set_xticklabels(plotorder)
        ax.set_title(title)
    else:
        plt.ylim(6,11) ###
        plt.boxplot(plotlist)
        plt.xticks(idxs,plotorder)
        plt.title(title)

    p = plt if not ax else ax
    for i in range(len(plotlist)):
        y = plotlist[i]
        x = np.random.normal(i+1, 0.04, size=len(y))
        p.plot(x, y, 'r.', markersize=10, alpha=0.5, c=colororder[i])

    if not ax:
        plt.savefig(path)
        plt.clf() # clear canvas

def plot_ori_inconsistency(indivsum_df, twosites_df, lbldf=False, namecol="Name", oricol='ori',
                affcol="affinity", log=False, fixed_ax=False, prefix_path = "", thresline=False):
    """
    Args:
        indivsum_df: data frame of Name, orientation, affinity of the individual sites sum
        twosites_df: data frame of Name, orientation, affinity of the two sites
        lbldf: data frame of Name, label in first orientation, label in second orientation
        cutoffline: float, plot cutoff line if specified
    """
    o1, o2 = indivsum_df[oricol].unique() # we assume we can only have 2 ori here
    indivsum, twosites = indivsum_df.copy(), twosites_df.copy()
    if log:
        indivsum[affcol] = np.log(indivsum[affcol])
        twosites[affcol] = np.log(twosites[affcol])
    indivsum["type"] = "one_" + indivsum[oricol]
    twosites["type"] = "two_" + twosites[oricol]

    numcol = 4
    numrow = 4
    if lbldf:
        label_names = set(lbldf[o1]) | set(lbldf[o2])
        perms = [(x,y) for x in label_names for y in label_names]
    else:
        perms = [["NA","NA"]]

    for perm in perms:
        if lbldf:
            cur_lbls = lbldf[(lbldf[o1] == perm[0]) & (lbldf[o2] == perm[1])]
            if cur_lbls.empty:
                continue
            cur_indiv = indivsum.merge(cur_lbls[[namecol]], on=namecol)
            cur_two = twosites.merge(cur_lbls[[namecol]], on=namecol)
        else:
            cur_indiv = indivsum
            cur_two = twosites
        curdf = pd.concat([cur_indiv, cur_two])
        print(perm, curdf[namecol].unique().shape[0])
        if len(curdf[affcol].tolist()) == 0:
            continue
        min_y = min(curdf[affcol].tolist())
        max_y = max(curdf[affcol].tolist())

        fig, ax = plt.subplots(numrow, numcol, figsize=(25, 5))
        plt.subplots_adjust(hspace = 0.4, wspace=0.2)
        all_lbl = "%s%s_%s_%s_%s" % (prefix_path,perm[0],o1,perm[1],o2)
        po = ["one_%s"%o1,"two_%s"%o1,"one_%s"%o2,"two_%s"%o2]
        co = ["red","red","blue","blue"]
        with PdfPages("%s.pdf" % all_lbl) as pdf:
            fig = plt.figure(figsize=(25,14))
            fig.subplots_adjust(hspace=0.4,wspace=0.2)
            i = 1 # start from 1 for the all plot
            cur_ax = fig.add_subplot(numcol,numrow,i)
            df_median = curdf.groupby([namecol,"type",oricol]).median().reset_index()
            scatter_boxplot_col(df_median, "type", affcol , plotorder=po, colororder=co, title=all_lbl, ax=cur_ax)
            for n in curdf[namecol].drop_duplicates().tolist():
                if i == 0:
                    fig = plt.figure(figsize=(25,14))
                    fig.subplots_adjust(hspace=0.4,wspace=0.2)
                i += 1
                cur_ax = fig.add_subplot(numcol,numrow,i)
                curdf_col = curdf[curdf[namecol] == n]
                scatter_boxplot_col(curdf_col, "type", affcol, plotorder=po, colororder=co, title=n, ax=cur_ax)
                # statistical annotation
                for pos in [1, 3]:
                    p = st.wilcox(curdf_col[curdf_col["type"] == po[pos-1]]["affinity"].tolist(),
                                  curdf_col[curdf_col["type"] == po[pos]]["affinity"].tolist(), 'less')
                    x1, x2 = pos, pos+1
                    y, h, col = max(curdf_col["affinity"]) + 0.2, 0.2, 'k'
                    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], c=col)
                    plt.text((x1+x2)*.5, y+h, "%.4f" % p, ha='center', va='bottom', color=col)
                if thresline:
                    cur_ax.axhline(thresline, linestyle='--', color="orange", lw=0.5)
                if fixed_ax:
                    cur_ax.set_ylim([min_y, max_y])
                if i == numcol*numrow:
                    pdf.savefig(fig)
                    plt.close()
                    i = 0
            pdf.savefig(fig)
    plt.clf()
    plt.close()
    return 0

#def plot_label_change()
