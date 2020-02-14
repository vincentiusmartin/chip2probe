'''
Created on Oct 30, 2019

@author: vincentiusmartin
'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from matplotlib.backends.backend_pdf import PdfPages
from decimal import Decimal

import trainingdata.seqextractor as seqextractor

import math

import statistics
import util.stats as st
import util.bio as bio

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
            self.df = trainingdata.copy().reset_index(drop=True)
        elif isinstance(trainingdata, str):
            self.df = pd.read_csv(trainingpath,sep=sep)
        else:
            raise Exception("input must be string or data frame")
        self.motiflen = corelen

    def get_labels_indexes(self):
          return self.df.groupby("label").groups

    def stacked_bar_categories(self, x, y=["label"],plotname="stackedbar.png",avg=False,legend=True, ylabel=""):
        cat_df = self.df[[x]]
        cat_df["label"] = self.df['label']
        group = [x] + y
        df2 = cat_df.groupby(group)['label'].count() # .unstack(x).fillna(0)
        if avg:
            df2 = df2.groupby(level=0).apply(lambda x: x / float(x.sum()))
        ax = df2.unstack(x).fillna(0).T.plot(kind='bar', stacked=True,legend=legend,rot=0)
        ax.set_ylabel(ylabel)
        plt.savefig(plotname)
        plt.clf()

    # SHOULDN'T BE HERE, TODO: MOVE SOMEWHERE ELSE
    def boxplot_categories(self, df, by=["label"], input_cols="default", plotname="boxplot.png"):
        if input_cols == "default":
            cols = list(set(df.columns) - set(by))
        else:
            cols = list(input_cols)
        sns.set_style("whitegrid")
        numcol = 4
        numrow = math.ceil(len(cols) / numcol)
        # to make axis with different y-scale
        fig, ax = plt.subplots(numrow, numcol, figsize=(14, 5))
        plt.subplots_adjust(hspace = 0.4, wspace=0.6)
        grouped = df.groupby(by=by)
        # need to sort to keep the order consistent
        cols.sort()
        for i in range(len(cols)):
            colname = cols[i]
            cur_group = {elm[0]:list(elm[1]) for elm in grouped[colname]}
            #cur_group = {elm[0]:list(filter(lambda a: a != 0, list(elm[1]))) for elm in grouped[colname]}
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
                p_gr = st.wilcox(data[x1],data[x2],alternative="two.sided")
                hfrac = -1 if mval < 0 else 1
                hfactor = (max2 - max1)**1.5
                pline_h = mval * 0.05
                pline_pos = mval * 0.2
                y, h, col = df[colname].max() + hfactor * pline_pos, pline_h, 'k'
                #cur_ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
                print(colname,"add>coop pvalue: ",p_gr)
                pstr = "%.2E" % Decimal(p_gr) if p_gr < 0.001 else "%.4f" % p_gr
                #cur_ax.text((x1+x2)*.5, y + h, "p = %s"%(pstr), ha='center', va='bottom', color="red")

                adj_ylim = (y+h) * 1.1
                if adj_ylim > ylim: ylim = adj_ylim
            #cur_ax.set_ylim(top=ylim)
        for d in range(len(cols),numrow*numcol):
            fig.delaxes(ax.flatten()[d])
        plt.savefig(plotname)
        plt.clf()

    def get_ordered_site_list(self):
        site1 = {}
        site2 = {}
        for idx,row in self.df.iterrows():
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1[idx] = row["site_str_pos"]
                site2[idx] = row["site_wk_pos"]
            else:
                site1[idx] = row["site_wk_pos"]
                site2[idx] = row["site_str_pos"]
        return site1,site2

    def training_summary(self, by=["label"], cols="default", plotname="training_summary.png"):
        if cols == "default":
            col_to_box = list(set(self.df.columns) - {"id", "name", "sequence", "label", "index"})
        else:
            col_to_box = cols
        self.boxplot_categories(self.df, by=["label"], input_cols=col_to_box, plotname="training_summary.png")

    def flip_one_face_orientation(self,positive_cores):
        # flip if orientation one
        ori = self.get_feature_orientation(positive_cores)
        records = self.df.to_dict(orient='records')
        for i in range(0,len(records)):
            if ori[i]["ori"] == "1":
                site_str = records[i]["sequence"][records[i]["site_str_pos"] - self.motiflen//2:records[i]["site_str_pos"] + self.motiflen//2]
                site_wk = records[i]["sequence"][records[i]["site_wk_pos"] - self.motiflen//2:records[i]["site_wk_pos"] + self.motiflen//2]
                if site_str not in positive_cores and site_wk not in positive_cores:
                    records[i]["sequence"] = bio.revcompstr(records[i]["sequence"])
                    # flip the position as well
                    str_pos = records[i]["site_str_pos"]
                    wk_pos = records[i]["site_wk_pos"]
                    records[i]["site_str_pos"] = len(records[i]["sequence"]) - wk_pos
                    records[i]["site_wk_pos"] = len(records[i]["sequence"]) - str_pos
                    #site_str_new = records[i]["sequence"][records[i]["site_str_pos"] - self.motiflen//2:records[i]["site_str_pos"] + self.motiflen//2]
                    #site_wk_new = records[i]["sequence"][records[i]["site_wk_pos"] - self.motiflen//2:records[i]["site_wk_pos"] + self.motiflen//2]
                    #if site_str_new not in positive_cores or site_wk_new not in positive_cores:
                    #    print(records[i]["sequence"],records[i]["site_str_pos"],records[i]["site_wk_pos"])
                    #    print(site_str, site_wk, site_str_new, site_wk_new)
                    #    break
        return Training(pd.DataFrame(records), corelen=self.motiflen)


    def weak_type_to_int(self,name):
        if name.endswith("_weak_s1"):
            return 1
        elif name.endswith("_weak_s2"):
            return 2
        elif name.endswith("_weak_s1_2"):
            return 3
        return 0

    def recur_dictify(self, frame):
        # https://stackoverflow.com/questions/19798112/convert-pandas-dataframe-to-a-nested-dict
        if len(frame.columns) == 1:
            if frame.values.size == 1: return frame.values[0][0]
            return frame.values.squeeze()
        grouped = frame.groupby(frame.columns[0])
        d = {k: self.recur_dictify(g.ix[:,1:]) for k,g in grouped}
        return d

    def plot_grouped_label(self, df, by, column, figname="boxplot.png"):
        g = self.recur_dictify(df[by + [column]])
        for lbl in g:
            l1 = g[lbl]["additive"].tolist()
            l2 = g[lbl]["cooperative"].tolist()
            p = st.t_test(l1,l2,alternative="greater")
            pstr = "%.2E" % Decimal(p) if p < 0.001 else "%.4f" % p
            print("ori %s additive > cooperative, p: %s" % (lbl,pstr))
        df.boxplot(by=by, column=column)
        plt.savefig(figname)
        plt.close()


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
                """
                if group["id"].str.endswith(r"o1|o2"):
                    g_o1 = group[group["id"].str.endswith("o1")]
                    g_o2 = group[group["id"].str.endswith("o2")]
                    curgroup = g_o1 if len(g_o1) > len(g_o2) else g_o2
                """
                curgroup = group
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
        # change this to -1 or -2 depend # HARDCODE NOW
        #dcopy["Group"] = dcopy.apply(lambda row: "_".join(row["name"].split("_")[:-2]), axis=1)
        dcopy["Group"] = dcopy.apply(lambda row: "_".join(row["name"].split("_")[:-1]), axis=1)
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
                """
                if group["id"].str.endswith(r"o1|o2"):
                    g_o1 = group[group["id"].str.endswith("o1")]
                    g_o2 = group[group["id"].str.endswith("o2")]
                    curgroup = g_o1 if len(g_o1) > len(g_o2) else g_o2
                """
                curgroup = group
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

    def extract_positional(self,seq, maxk = 2, label="seq", minseqlen=-float("inf")):
        '''
        orientation: if right, then start from 0 to the right, else start from
        len(seq)-1 to the left
        minseqlen: will fill -1 if not enough bases
        '''
        iterseq = str(seq)
        nucleotides = ['A','C','G','T']
        features = {}
        for k in range(1,maxk+1):
            perm = ["".join(p) for p in itertools.product(nucleotides, repeat=k)]
            i = 0
            while i < len(iterseq)+1-k:
                for kmer in perm:
                    seqcmp = iterseq[i:i+k]
                    if seqcmp == kmer:
                        features["%s_pos%d_%s" % (label,i,kmer)] = 1
                    else:
                        features["%s_pos%d_%s" % (label,i,kmer)] = 0
                i += 1
            # append the rest with -1
            if minseqlen > 0:
                while i < minseqlen + 1 - k:
                    for kmer in perm:
                        features["%s_pos%d_%s" % (label,i,kmer)] = -1
                    i += 1
        return features
    # =========

    def get_feature_site_pref(self):
        rfeature = []
        for idx,row in self.df.iterrows():
            f = {"site_wk_score":row["site_wk_score"], "site_str_score":row["site_str_score"]}
            rfeature.append(f)
        return rfeature

    def get_nonrev_dfeature(self, k, label="seq", initcount=0):
        nucleotides = ['A','C','G','T']
        perm = []
        dfeature = {}
        for p in itertools.product(nucleotides, repeat=k):
            p = "".join(p)
            rev = bio.revcompstr(p)
            if p not in perm and rev not in perm:
                to_append = p if p < rev else rev
                dfeature["%s_%s"%(label,to_append)] = initcount
                perm.append(p)
        return dfeature

    def count_nonrev_kmer(self,seq,maxk,label="seq",avg=False):
        dfeature = {}
        for k in range(1, maxk+1):
            dfeature = {**dfeature,**self.get_nonrev_dfeature(k,label,initcount=0)}
            for i in range(len(seq) + 1 - k):
                curseq = str(seq)[i:i+k]
                revseq = bio.revcompstr(curseq)
                curseq = curseq if curseq < revseq else revseq
                key = "%s_%s"%(label,curseq)
                dfeature[key] = dfeature[key] + 1
        if avg:
            for key in dfeature:
                kmer = key[len("%s_"%label):]
                if len(seq) + 1 - len(kmer) != 0:
                    dfeature[key] = dfeature[key] / (len(seq) + 1 - len(kmer))
        return dfeature

    def get_middle_avgshape_feature(self, freqs, dnashape, maxk=2, action="avg"):
        rfeature = []
        for idx,row in self.df.iterrows():
            dfeature = {}
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
            # since position is the middle point of each site
            start = site1 + self.motiflen // 2
            end = site2 - self.motiflen // 2
            shapes = {"prot":dnashape.prot, "mgw":dnashape.mgw, "roll":dnashape.roll, "helt":dnashape.helt}
            #shapes = {"mgw":dnashape.mgw}
            linker_len = end - start
            for s in shapes:
                linker_shape = shapes[s][str(idx + 1)][start:end]
                for freq in freqs:
                    span = freq // 2
                    label = "middle_%s_%s_%d_%d" % (s, action, freq, freq + 1)
                    if linker_len < freq:
                        for k in range(0,maxk):
                            midseqval = [-999]
                    else:
                        mid = linker_len // 2
                        if linker_len % 2 == 0: # even
                            midseqval = linker_shape[mid-span-1:mid+span+1]
                        else: # odd
                            midseqval = linker_shape[mid-span:mid+span+1]
                    if action == "mean":
                        midseqval = statistics.mean(midseqval)
                    elif action == "max":
                        midseqval = max(midseqval)
                    elif action == "min":
                        midseqval = min(midseqval)
                    dfeature[label] = midseqval
            rfeature.append(dfeature)
        return rfeature

    def get_middle_feature(self, freqs, maxk=2):
        # pos needs to be odd number
        # the representation is average
        rfeature = []
        for idx,row in self.df.iterrows():
            dfeature = {}
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
            # since position is the middle point of each site
            start = site1 + self.motiflen // 2
            end = site2 - self.motiflen // 2
            linker = row["sequence"][start:end]
            for freq in freqs:
                span = freq // 2
                label = "middle_%d_%d" % (freq, freq + 1)
                if len(linker) < freq:
                    for k in range(0,maxk):
                        dfeature = {**dfeature,**self.get_nonrev_dfeature(k + 1, label=label, initcount=-1)}
                else:
                    mid = len(linker) // 2
                    if len(linker) % 2 == 0: # even
                        midseq = linker[mid-span-1:mid+span+1]
                    else: # odd
                        midseq = linker[mid-span:mid+span+1]
                    dfeature = {**dfeature, **self.count_nonrev_kmer(midseq,maxk=maxk,label=label,avg=True)}
            rfeature.append(dfeature)
        return rfeature

    def get_feature_flank_shapes(self, dnashape, seqin):
        # seqin: how many base inside the linker
        rfeature = []
        for idx,row in self.df.iterrows():
            dfeature = {}
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
                s1type, s2type = "str", "wk"
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
                s1type, s2type = "wk", "str"
            # since position is the middle point of each site
            start = site1 + self.motiflen // 2
            end = site2 - self.motiflen // 2
            shapes = {"prot":dnashape.prot, "mgw":dnashape.mgw, "roll":dnashape.roll, "helt":dnashape.helt}
            #shapes = {"mgw":dnashape.mgw}
            linker_len = end - start
            minlink = seqin if seqin < linker_len else linker_len
            for s in shapes:
                linker_shape = shapes[s][str(idx + 1)][start:end]
                linker1 = linker_shape[:minlink]
                linker2 = linker_shape[-minlink:][::-1]
                i = 0
                while i < minlink:
                    dfeature["%s_%s_pos-%d" % (s,s1type,i)] = linker1[i]
                    dfeature["%s_%s_pos-%d" % (s,s2type,i)] = linker2[i]
                    i += 1
                while i < seqin:
                    dfeature["%s_%s_pos-%d" % (s,s1type,i)] = -999
                    dfeature["%s_%s_pos-%d" % (s,s2type,i)] = -999 # currently 0 ???
                    i += 1
            rfeature.append(dfeature)
        return rfeature

    # only linker now
    def get_feature_flank_seqs(self, k, seqin=0):
        rfeature = []
        for idx,row in self.df.iterrows():
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
                s1type, s2type = "str", "wk"
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
                s1type, s2type = "wk", "str"
            # since position is the middle point of each site
            start = site1 + self.motiflen // 2
            end = site2 - self.motiflen // 2
            linker = row["sequence"][start:end]
            minlink = seqin if seqin < len(linker) else len(linker)
            flank1 = linker[:minlink]
            flank2 = linker[-minlink:][::-1] # also reverse the string
            d1 = self.extract_positional(flank1,maxk=k,minseqlen=seqin,label="flankseq_%s"%s1type)
            d2 = self.extract_positional(flank2,maxk=k,minseqlen=seqin,label="flankseq_%s"%s2type)
            rfeature.append({**d1, **d2})
        return rfeature

    def get_feature_distance(self, type="numerical"):
        if type == "numerical":
            return [{"dist-numeric":x} for x in self.df["distance"].values]
        elif type == "categorical":
            one_hot = pd.get_dummies(self.df['distance'])
            one_hot.columns = ["dist-cat-%d"%col for col in one_hot.columns]
            return one_hot.to_dict('records')
        else:
            raise Exception("distance must be numerical or categorical")

    def get_linker_list(self):
        linkers = []
        falses = []
        for idx,row in self.df.iterrows():
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
            # since position is the middle point of each site
            start = site1 + self.motiflen // 2
            end = site2 - self.motiflen // 2
            if site1 + site2 != 36 and site1 + site2 != 37:
                falses.append(row["sequence"])
            linker = row["sequence"][start:end]
            linkers.append({"linker":linker, "site1":site1, "site2":site2})
        #print(len(falses))
        return linkers

    def get_linker_GC_content(self):
        rfeature = []
        for idx,row in self.df.iterrows():
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
            # since position is the middle point of each site
            start = site1 + self.motiflen // 2
            end = site2 - self.motiflen // 2
            linker = row["sequence"][start:end]
            if len(linker) == 0:
                rfeature.append({"GC_content":0})
            else:
                rfeature.append({"GC_content":float(linker.count('G') + linker.count('C'))/len(linker)})
        return rfeature

    def get_feature_linker_composition(self, k):
        rfeature = []
        for idx,row in self.df.iterrows():
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
            # since position is the middle point of each site
            start = site1 + self.motiflen // 2
            end = site2 - self.motiflen // 2
            linker = row["sequence"][start:end]
            ratio = seqextractor.extract_kmer_ratio(linker,k)
            rfeature.append(ratio)
        return rfeature

    def get_feature_orientation(self, positive_cores, relative=True, one_hot=False):
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
                    ori = 'HT/TH'
                elif s1 == -1 and s2 == -1:
                    ori = 'HT/TH'
                elif s1 == 1 and s2 == -1:
                    ori = 'HH'
                elif s1 == -1 and s2 == 1:
                    ori = 'TT'
                else:
                    ori = '-1'
                rfeature.append({"ori":ori})
            else:
                rfeature.append({"ori1":s1, "ori2":s2})
        if one_hot:
            return pd.get_dummies(pd.DataFrame(rfeature)).to_dict('records')
        else:
            return rfeature
