import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools

import chip2probe.util.stats_r as st

def make_replicas_permutation(df, oricol="ori", namecol="Name", affcol="affinity", typecol="type", repcol="rep"):
    """
    Make permutation across replicas
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
            print(g)
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
            break
        break
    return indivsum, twosites

def create_cooplbl(twosites, indivsum, pcutoff = 0.05):
    """
    """
    p_coop = st.wilcox(twosites, indivsum, "greater")
    p_anti = st.wilcox(twosites, indivsum, "less")
    if p_coop < pcutoff:
        return "cooperative", p_coop
    elif p_coop < pcutoff:
        return "anticoop", p_coop
    else:
        return 'additive', p_coop

# --------

def plot_chamber_corr(dfx, dfy, xlab="Chamber 1", ylab="Chamber 2",
                namecol="Name", valcol="Alexa488Adjusted", repcol="rep", extrajoincols=[],
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
        Show plot result for correlation between chamber
    """
    # correlation plot of negative controls in chamber 1 vs chamber 2

    joincols = [namecol, repcol] + extrajoincols
    dfcombined = dfx.merge(dfy, on=joincols)[joincols + ["%s_x"%valcol, "%s_y"%valcol]]

    if median:
        dfcombined = dfcombined.groupby([namecol] + extrajoincols)[["%s_x"%valcol, "%s_y"%valcol]].median().reset_index()
    if log:
        dfcombined["%s_x"%valcol] = np.log(dfcombined["%s_x"%valcol])
        dfcombined["%s_y"%valcol] = np.log(dfcombined["%s_y"%valcol])
        cutoff = np.log(cutoff)

    x = dfcombined["%s_x"%valcol].values
    y = dfcombined["%s_y"%valcol].values

    slope, intercept = np.polyfit(x, y, 1)
    # Create a list of values in the best fit line
    abline_values = [slope * i + intercept for i in x]

    f = plt.figure()
    ax = f.add_subplot(111)
    # plot diagonal
    plt.plot([min(min(x),min(y)),max(max(x),max(y))], [min(min(x),min(y)),max(max(x),max(y))], color='black', label='diagonal', linewidth=0.5)

    if cutoff:
        plt.axhline(cutoff, color='black', linewidth=0.5, linestyle=":")
        plt.axvline(cutoff, color='black', linewidth=0.5, linestyle=":")

    # plot best fit line
    plt.plot(x, abline_values, color='red', label='best fit', linewidth=0.5)
    r_squared = np.corrcoef(x,y)[0,1]**2
    plt.scatter(x, y, s=3)

    if shownames:
        names = dfcombined[namecol].values
        for i in ranplot_classified_labelsge(len(names)):
            plt.text(x[i], y[i], names[i], fontsize=3)

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.text(0.02, 0.94, 'R² = %0.2f' % r_squared, color='red', transform = ax.transAxes)
    sign = '+' if intercept > 0 else '-'
    plt.text(0.02, 0.9, 'best fit: y=%0.2fx %s %0.2f' % (slope,sign,abs(intercept)), color='red', transform = ax.transAxes)
    plt.legend(loc='lower right')

    if path:
        plt.savefig(path)
    else:
        plt.show()
    plt.clf()

def plot_classified_labels(df, path="", col1="Alexa488Adjusted_x", col2="Alexa488Adjusted_y",
                           xlab="Chamber1", ylab="Chamber2", log=True, title=""):
    """
    Desc

    Args:

    Return:

    """
    lblclr = [("below_cutoff","yellow", 0.7), ("additive","gray",1), ("cooperative","blue",1),("anticoop","red",1)]
    newdf = df.copy()
    if log:
        newdf[col1] = np.log(df[col1])
        newdf[col2] = np.log(df[col2])

    for lc in lblclr:
        if not newdf[newdf['label']==lc[0]].empty:
            x = newdf[newdf['label']==lc[0]][col1].values
            y = newdf[newdf['label']==lc[0]][col2].values
            label = lc[0] if lc[0] != "below_cutoff" else None
            plt.scatter(x, y, color=lc[1], s=3, alpha=lc[2], label=label)

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    x,y = newdf[col1].values, newdf[col2].values
    plt.plot([min(min(x),min(y)),max(max(x),max(y))], [min(min(x),min(y)),max(max(x),max(y))], color='black')
    plt.legend(loc="lower right")

    plt.title(title)
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

def scatter_boxplot_col(input_df, coltype, colval, plotorder=[], ax=None, path = "box.png", title=""):
    plotlist = []
    for po in plotorder:
        plotlist.append(input_df[input_df[coltype] == po][colval].tolist())

    idxs = list(range(1,len(plotorder)+1))
    if ax:
        ax.boxplot(plotlist)
        ax.set_xticklabels(plotorder)
        ax.set_title(title)
    else:
        plt.boxplot(plotlist)
        plt.xticks(idxs,plotorder)
        plt.title(title)

    p = plt if not ax else ax
    for i in range(len(plotlist)):
        y = plotlist[i]
        x = np.random.normal(i+1, 0.04, size=len(y))
        p.plot(x, y, 'r.', alpha=0.2)

    if not ax:
        plt.savefig(path)
        plt.clf() # clear canvas

def plot_ori_incosistency2(indivsum_df, twosites_df, lbldf, namecol="Name",
                affcol="affinity", log=False, fixed_ax=False):
    """
    labels: data frame
    """
    o1, o2 = set(lbldf.columns) - {namecol}
    indivsum, twosites = indivsum_df.copy(), twosites_df.copy()
    if log:
        indivsum[affcol] = np.log(indivsum[affcol])
        twosites[affcol] = np.log(twosites[affcol])
    indivsum["type"] = "one_" + indivsum["ori"]
    twosites["type"] = "two_" + twosites["ori"]

    numcol = 4
    numrow = 4
    label_names =  ["cooperative", "additive", "anticoop"]
    perms = [(x,y) for x in label_names for y in label_names]

    for perm in perms:
        cur_lbls = lbldf[(lbldf[o1] == perm[0]) & (lbldf[o2] == perm[1])]
        if cur_lbls.empty:
            continue
        cur_indiv = indivsum.merge(cur_lbls[["Name"]], on="Name")
        cur_two = twosites.merge(cur_lbls[["Name"]], on="Name")
        curdf = pd.concat([cur_indiv, cur_two])
        print(perm, curdf["Name"].unique().shape[0])
        min_y = min(curdf[affcol].tolist())
        max_y = max(curdf[affcol].tolist())

        fig, ax = plt.subplots(numrow, numcol, figsize=(25, 5))
        plt.subplots_adjust(hspace = 0.4, wspace=0.2)
        all_lbl = "%s_%s_%s_%s" % (perm[0],o1,perm[1],o2)
        po = ["one_%s"%o1,"two_%s"%o1,"one_%s"%o2,"two_%s"%o2]
        with PdfPages("%s.pdf" % all_lbl) as pdf:
            fig = plt.figure(figsize=(25,14))
            fig.subplots_adjust(hspace=0.4,wspace=0.2)
            i = 1 # start from 1 for the all plot
            cur_ax = fig.add_subplot(numcol,numrow,i)
            df_median = curdf.groupby(["Name","type","ori"]).median().reset_index()
            scatter_boxplot_col(df_median, "type", "affinity" , plotorder=po, title=all_lbl, ax=cur_ax)
            for n in curdf["Name"].drop_duplicates().tolist():
                if i == 0:
                    fig = plt.figure(figsize=(25,14))
                    fig.subplots_adjust(hspace=0.4,wspace=0.2)
                i += 1
                cur_ax = fig.add_subplot(numcol,numrow,i)
                scatter_boxplot_col(curdf[curdf["Name"] == n], "type", "affinity", plotorder=po, title=n, ax=cur_ax)
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

def plot_ori_inconsistency(df1, df2, cols, orinames=["o1","o2"], log=False, fixed_ax=True):
    """
    Desc

    Args:
        fixed_ax: have to lower and upper y axis value to be the same for all plots
    Return:

    """
    o1,o2 = orinames
    df = df1.merge(df2, on=["Name", "rep"], suffixes=('_%s'%o1, '_%s'%o2))
    oricols = ["%s_%s" % (c,o) for o in orinames for c in cols ]
    if log:
        for c in oricols:
            df[c] = np.log(df[c])

    numcol = 4
    numrow = 4
    labels =  ["cooperative", "additive", "anticoop"]
    perms = [(x,y) for x in labels for y in labels]

    for perm in perms:
        curdf = df[(df["label_%s"%o1] == perm[0]) & (df["label_%s"%o2] == perm[1])]
        print(perm, curdf["Name"].unique().shape[0])
        min_y, max_y = min(df.loc[:,oricols].min().tolist()), max(df.loc[:,oricols].max().tolist())
        all_lbl = "%s_%s_%s_%s" % (perm[0],o1,perm[1],o2)
        fig, ax = plt.subplots(numrow, numcol, figsize=(25, 5))
        plt.subplots_adjust(hspace = 0.4, wspace=0.2)
        with PdfPages("%s.pdf" % all_lbl) as pdf:
            fig = plt.figure(figsize=(25,14))
            fig.subplots_adjust(hspace=0.4,wspace=0.2)
            i = 1 # start from 1 for the all plot
            cur_ax = fig.add_subplot(numcol,numrow,i)
            dfall_median = curdf.groupby("Name").median().reset_index()
            print(dfall_median)
            scatter_boxplot(dfall_median, cols=oricols, title=all_lbl, ax=cur_ax)
            for n in curdf["Name"].drop_duplicates().tolist():
                if i == 0:
                    fig = plt.figure(figsize=(25,14))
                    fig.subplots_adjust(hspace=0.4,wspace=0.2)
                i += 1
                cur_ax = fig.add_subplot(numcol,numrow,i)
                scatter_boxplot(curdf[curdf["Name"] == n], cols=oricols, log=True, title=n, ax=cur_ax)
                if fixed_ax:
                    cur_ax.set_ylim([min_y, max_y])
                if i == numcol*numrow:
                    pdf.savefig(fig)
                    plt.close()
                    i = 0
            pdf.savefig(fig)
        plt.clf()
        plt.close()
