'''
Created on Oct 30, 2019

@author: vincentiusmartin
'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from matplotlib.backends.backend_pdf import PdfPages

import trainingdata.seqextractor as seqextractor

import math

from statistics import median
import util.stats as st

class Training(object):
    '''
    classdocs
    '''

    # TODO: make column name more general
    def __init__(self, trainingdata, corelen, sep="\t"):
        '''
        Constructor
        '''
        if isinstance(trainingdata, pd.DataFrame):
            self.df = trainingdata.copy()
        elif isinstance(trainingdata, str):
            self.df = pd.read_csv(trainingpath,sep=sep)
        else:
            raise Exception("input must be string or data frame")
        self.motiflen = corelen


    def plot_summary(self, by=["label"], cols="default", plotname="training_summary.png"):
        if cols == "default":
            col_to_box = list(set(self.df.columns) - {"id", "name", "sequence", "label"})
        else:
            col_to_box = cols
        sns.set_style("whitegrid")
        numcol = 3
        numrow = math.ceil(len(col_to_box) / numcol)
        # to make axis with different y-scale
        fig, ax = plt.subplots(numrow, numcol, figsize=(14, 5))
        plt.subplots_adjust(hspace = 0.4, wspace=0.6)
        grouped = self.df.groupby(by=by)
        # need to sort to keep the order consistent
        col_to_box.sort()
        for i in range(len(col_to_box)):
            colname = col_to_box[i]
            cur_group = {elm[0]:list(elm[1]) for elm in grouped[colname]}
            labels, data = [*zip(*cur_group.items())]
            #labels = list(cur_group.keys())#["cooperative","additive","anticoop"]
            #data = [cur_group[x] for x in labels]
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
                p_gr = st.wilcox(data[x1],data[x2],alternative="greater")
                hfactor = (x2 - x1)**2.1
                pline_h = mval * 0.1
                pline_pos = mval * 0.05
                y, h, col = self.df[colname].max() + hfactor * pline_pos, pline_h, 'k'
                cur_ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
                cur_ax.text((x1+x2)*.5, y + h, "p-gr=%.3f"%(p_gr), ha='center', va='bottom', color="red")

                adj_ylim = (y+h) * 1.1
                if adj_ylim > ylim: ylim = adj_ylim
            cur_ax.set_ylim(top=ylim)
        for d in range(len(col_to_box),numrow*numcol):
            fig.delaxes(ax.flatten()[d])
        plt.savefig(plotname)
        plt.clf()

    def weak_type_to_int(self,name):
        if name.endswith("_weak_s1"):
            return 1
        elif name.endswith("_weak_s2"):
            return 2
        elif name.endswith("_weak_s1_2"):
            return 3
        return 0

    def plot_weak_sites(self, filepath="weak.pdf"):
        dcopy = self.df.copy()
        dcopy["Group"] = dcopy.apply(lambda row: row["name"].split("_weak")[0], axis=1)
        dcopy["label"] = dcopy['label'].map({'cooperative': 1, 'additive': 0, 'anticoop':-1})
        dcopy["type"] = dcopy.apply(lambda row: self.weak_type_to_int(row["name"]), axis=1)
        #df_o1 = dcopy[dcopy["id"].str.endswith("o1")]
        groups = dcopy.groupby("Group")
        numcol = 4
        numrow = 4
        fig, ax = plt.subplots(numrow, numcol, figsize=(25, 5))
        plt.subplots_adjust(hspace = 0.4, wspace=0.2)
        with PdfPages(filepath) as pdf:
            fig = plt.figure(figsize=(25,14))
            fig.subplots_adjust(hspace=0.4,wspace=0.4)
            n = 0
            for name, group in groups:
                g_o1 = group[group["id"].str.endswith("o1")]
                g_o2 = group[group["id"].str.endswith("o2")]
                curgroup = g_o1 if len(g_o1) > len(g_o2) else g_o2
                if len(curgroup) > 2:
                    if n == 0:
                        fig = plt.figure(figsize=(25,14))
                        fig.subplots_adjust(hspace=0.4,wspace=0.2)
                    n += 1
                    ax = fig.add_subplot(numcol,numrow,n)
                    ax = sns.scatterplot(x='type', y='label', data=curgroup)
                    ax.set_xticks([0,1,2,3])
                    ax.set_xticklabels(['wt','weak_s1','weak_s2','weak_both']) # set the labels
                    ax.set_yticks([-1,0,1])
                    ax.set_yticklabels(['anticoop','additive','cooperative']) # set the labels
                    ax.set_ylabel("")
                    ax.set_ylim(-2,2) # so we have everything in the middle
                    ax.set_title(name)
                    if n == numcol*numrow:
                        pdf.savefig(fig)
                        plt.close()
                        n = 0
            pdf.savefig(fig)
        plt.clf()
        plt.close()

    """
    def plot_by_group():
        g = self.df.groupby('site%d_pos' % site)
        grouped = {str(key):item for key,item in g}
        sns.countplot(x='distance', y='label', data=self.df)
        plt.savefig(plotname)
        plt.clf()
    """

    # right now, this is only for custom distance
    def plot_distance_numeric(self, filepath="dist.pdf"):
        dcopy = self.df.copy()
        dcopy["Group"] = dcopy.apply(lambda row: "_".join(row["name"].split("_")[:-2]), axis=1)
        dcopy["label"] = dcopy['label'].map({'cooperative': 1, 'additive': 0, 'anticoop':-1})
        #df_o1 = dcopy[dcopy["id"].str.endswith("o1")]
        groups = dcopy.groupby("Group")

        numcol = 4
        numrow = 4
        fig, ax = plt.subplots(numrow, numcol, figsize=(25, 5))
        plt.subplots_adjust(hspace = 0.4, wspace=0.2)
        with PdfPages(filepath) as pdf:
            fig = plt.figure(figsize=(25,14))
            fig.subplots_adjust(hspace=0.4,wspace=0.4)
            n = 0
            progress = 0
            modprog = len(groups) // 10 + 1
            for name, group in groups:
                if progress % modprog == 0:
                    print("Progress: %.1f%%" % (100 * float(progress)/len(groups)))
                progress += 1
                g_o1 = group[group["id"].str.endswith("o1")]
                g_o2 = group[group["id"].str.endswith("o2")]
                curgroup = g_o1 if len(g_o1) > len(g_o2) else g_o2
                if len(curgroup) > 2:
                    if n == 0:
                        fig = plt.figure(figsize=(25,14))
                        fig.subplots_adjust(hspace=0.4,wspace=0.2)
                    n += 1
                    ax = fig.add_subplot(numcol,numrow,n)
                    curgroup = curgroup.sort_values(by=['distance'])
                    curgroup["distance"] = curgroup["distance"].astype(str)
                    ax = sns.scatterplot(x='distance', y='label', data=curgroup)
                    ax.set_ylim(-2,2)
                    ax.set_yticks([-1,0,1])
                    ax.set_yticklabels(['anticoop','additive','cooperative']) # set the labels
                    ax.set_ylabel("")
                    ax.set_title(name)
                    if n == numcol*numrow:
                        pdf.savefig(fig)
                        plt.close()
                        n = 0
            print("Progress: %.1f%%" % (100 * float(progress)/len(groups)))
            pdf.savefig(fig)
        plt.clf()
        plt.close()
    # =========

    def get_feature_site_pref(self):
        rfeature = []
        for idx,row in self.df.iterrows():
            f = {"site_wk_score":row["site_wk_score"], "site_str_score":row["site_str_score"]}
            rfeature.append(f)
        return rfeature

    def get_feature_distance(self, type="numerical"):
        if type == "numerical":
            return [{"dist-numeric":x} for x in self.df["distance"].values]
        elif type == "categorical":
            one_hot = pd.get_dummies(self.df['distance'])
            one_hot.columns = ["dist-cat-%d"%col for col in one_hot.columns]
            return one_hot.to_dict('records')
        else:
            raise Exception("distance must be numeric or categorical")

    def get_feature_linker_composition(self, k):
        rfeature = []
        for idx,row in self.df.iterrows():
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
            # since position is the middle point of each site
            start = site1 + self.motiflen // 2 + 1
            end = site2 - self.motiflen // 2
            linker = row["sequence"][start:end]
            ratio = seqextractor.extract_kmer_ratio(linker,k)
            rfeature.append(ratio)
        return rfeature

    def get_feature_orientation(self, positive_cores, relative=True):
        negative_cores = [seqextractor.revcompstr(p) for p in positive_cores]
        rfeature = []
        for idx,row in self.df.iterrows():
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
            seq = row["sequence"]
            p1 = seq[site1 - self.motiflen//2:site1 + self.motiflen//2]
            p2 = seq[site2 - self.motiflen//2:site2 + self.motiflen//2]
            if p1 in positive_cores:
                s1 = 1
            elif p1 in negative_cores:
                s1 = -1
            else:
                s1 = 0
                print("couldn't find the first site %s in %s in the core list" % (p1,seq))

            if p2 in positive_cores:
                s2 = 1
            elif p2 in negative_cores:
                s2 = -1
            else:
                s2 = 0
                print("couldn't find the second site %s in %s in the core list" % (p2,seq))
            if relative:
                if s1 == 1 and s2 == 1:
                    ori = '0'
                elif s1 == -1 and s2 == -1:
                    ori = '1'
                elif s1 == 1 and s2 == -1:
                    ori = '2'
                elif s1 == -1 and s2 == 1:
                    ori = '3'
                else:
                    ori = '-1'
                rfeature.append({"ori":ori})
            else:
                rfeature.append({"ori1":s1, "ori2":s2})
        return rfeature
