
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import string, random
import pathlib
import operator

from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

from chip2probe.util import bio
from chip2probe.util import stats as st
from chip2probe.util import util

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

dnashape_r = importr('DNAshapeR')

class DNAShape:

    def __init__(self, fasta): #bpos1_list, bpos2_list
        """
        DNA Shape object

        Args:
            fasta: a dictionary of sequence label to sequence or a filepath to a
                fasta file

         Returns:
            NA
        """
        self.shapetypes = ["ProT" , "MGW", "Roll", "HelT"]
        if isinstance(fasta, dict) or isinstance(fasta, list):
            # need to make fasta file first because DNAShapeR only accepts filepath
            tmpstr = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
            path = "%s.fasta" % tmpstr
            bio.makefasta(fasta, path)
            ds = dnashape_r.getShape(path)
        elif isinstance(fasta, str):
            path = str(fasta)
        else:
            raise TypeError("fasta must be a dictionary or path to a fasta file")
        for st in self.shapetypes:
            shpath = "%s.%s"%(path,st)
            setattr(self, st.lower(), util.read_dictlist_file(shpath))

        if isinstance(fasta, dict) or isinstance(fasta, list):
            # remove the tmp files from DNAShapeR
            for p in pathlib.Path(".").glob("%s*" % tmpstr):
                p.unlink()

    def get_shape_names():
        return self.shapetypes

    def get_sequence_count(self, seqlist, avg=True):
        seqlen = len(seqlist[0])
        if not all(len(l) == seqlen for l in seqlist):
            raise ValueError('not all sequences have same length!')
        counter = []
        for i in range(0,seqlen):
            nuc_ct = {"A":0, "C":0, "G": 0, "T":0}
            for j in range(0, len(seqlist)):
                nuc_ct[seqlist[j][i]] += 1
            if avg:
                nuc_ct = {k: float(v)/len(seqlist) for k, v in nuc_ct.items()}
            counter.append(nuc_ct)
        return counter

    # hard coded for coop
    # assume: seq same length, same site position
    #def plot_average(self, labels, site1list, site2list, sequences, pthres=0.05, path="shape.pdf", plotlabel="Average DNA shape"):
    def plot_average(self, df, linemarks=[], path="shape.pdf", lblcol="label", seqcol="sequence",
             pthres=0.01, pltlabel=""):
        """
        Plot average DNA shape

        Args:
            labels: list
            sites_per_tf:
            pthres: p-value threshold for difference between label
         Returns:
            NA
        TODO:
            label is still hardcoded, make it general?
        """
        plt.style.use('seaborn-whitegrid')
        # when the sequence is not aligned, we get the most frequent s1 and do
        # alignment according to this

        shapes = {s:getattr(self,s.lower()) for s in self.shapetypes}

        # lblidxs.get(l, 0)

        # get the index of each label
        lblidxs = {}
        # get the count of the nucleotide for the mini subplot
        seqlen = len(df['sequence'].iloc[0])
        base_ct = {}
        for l in df[lblcol].unique():
            seqlist = df[df[lblcol]==l]["sequence"].tolist()
            base_ct[l] = self.get_sequence_count(seqlist, avg=False)

        fig = plt.figure(figsize=(12,12))
        labels = df[lblcol].unique()
        n = 0 # for the subplot
        colors = {"additive":'orangered',"cooperative":'dodgerblue'} # TODO: make this general

        keystolabel = {"additive":"ã","cooperative":"ç"}
        for sh in shapes:
            n += 1
            ax = fig.add_subplot(2,2,n)
            lshift_dict = {}
            flen = len(next(iter(shapes[sh].values())))
            xlist = list(range(seqlen-flen,seqlen))
            # ==== Plotting part using 25, 50, and 75 percentiles ====
            min_y = float("inf")
            yall = {}
            for label in labels:
                curnames = df[df["label"] == label]["Name"].tolist()
                shapes_lbl = operator.itemgetter(*curnames)(shapes[sh])
                if not all(len(l) == flen for l in shapes_lbl):
                    raise ValueError('not all lists have same length!')
                ylist = np.median(shapes_lbl, axis=0)
                y25p = np.percentile(shapes_lbl, 25, axis=0)
                y75p = np.percentile(shapes_lbl, 75, axis=0)
                ax.plot(xlist, ylist, alpha=0.8, label=label, c=colors[label], marker='o')
                ax.fill_between(xlist, y75p, y25p, alpha=0.15, facecolor=colors[label])
                yall[label] = shapes_lbl # save the y to do hypothesis testing later

            # ==== Mark binding sites as given from the input ====
            for m in linemarks:
                ax.axvline(x=m,linewidth=1, color='g',linestyle='--')

            # ==== Hypothesis testing to mark significant binding sites ====
            signiflabel = ['']*(seqlen-flen)
            if "cooperative" in labels and "additive" in labels:
                for i in range(0,flen):
                    # for now assume yall is of size 2
                    arr_coop = [y[i] for y in yall['cooperative']]
                    arr_add = [y[i] for y in yall['additive']]
                    if not any(np.isnan(x) for x in arr_coop) and not any(np.isnan(x) for x in arr_add):
                        p1 = st.wilcox(arr_coop,arr_add,"greater")
                        p2 = st.wilcox(arr_coop,arr_add,"less")
                        if p1 <= pthres:
                            signiflabel.append(keystolabel["cooperative"])
                        elif p2 <= pthres:
                            signiflabel.append(keystolabel["additive"])
                        else:
                            signiflabel.append('')
                    else:
                        signiflabel.append('')

            # #label = ['' if x not in signiflist else '*' for x in xi]
            ax.set_xticks(list(range(seqlen)))
            ax.set_xticklabels(signiflabel)
            ax.yaxis.set_label_text(sh)
            ax.xaxis.set_label_text('Sequence position')
            ax.legend(loc="upper right")
            ax.set_title(pltlabel if pltlabel else sh)
            low_y = ax.get_ylim()[0]
            hfactor = ax.get_ylim()[1] -  ax.get_ylim()[0]
            ax.set_ylim(bottom = low_y - np.abs(0.35 * hfactor))
            for i in range(0,4): # low_y here?
                ax.yaxis.get_major_ticks()[i].label1.set_visible(False)

            base_color = {'A': "red", 'C': "green" , 'G': "orange", 'T': "blue"}
            bot_anchor = 0
            for base in base_color:
                for key in base_ct:
                    inset_ax = inset_axes(ax,
                                  height="5%",
                                  width="100%",
                                  bbox_to_anchor= (0, bot_anchor, 1, 1),
                                  bbox_transform=ax.transAxes,
                                  loc=8)
                    inl = [elm[base] for elm in base_ct[key]]
                    # also add 0.5 here to center
                    xbar = [i+0.5 for i in range(0,len(inl))] # should use the same axis with ax but this works...
                    sns.barplot(x = xbar, y = inl, color = base_color[base], ax=inset_ax)
                    inset_ax.get_xaxis().set_visible(False)
                    # inset_ax.set_ylim(top=max)
                    inset_ax.set_yticklabels([])
                    inset_ax.patch.set_visible(False)
                    inset_ax.set_ylabel("%s_%s"%(base,key[:2]),rotation=0)
                    inset_ax.yaxis.label.set_color(base_color[base])
                    bot_anchor += 0.05
        with PdfPages(path) as pdf:
            pdf.savefig(fig)

# ============== This is a separate class to contain everything ==============

"""
class DNAShapes:

    def __init__(self, path, bsites):
        # INITIALIZE
        self.shapedict = {}
        self.dists = []
        self.bsites = bsites

        dirs = next(os.walk(path))[1]
        self.dists = [d[1:] for d in dirs]
        self.maxdist = max(self.dists)
        for distdir in dirs:
            distpath = "%s/%s" % (path,distdir)
            shape = DNAShape(distpath)
            self.shapedict[distdir] = shape

    def get_features(self, as_tbl=False):
        span_out = 1
        span_in = 5
        numseq = len(self.bsites[0])
        shape_names = DNAShape.get_shape_names()
        #features = {name:[[] for _ in range(numseq)] for name in shape_names}
        features = {}

        for dist,shape_dist in self.shapedict.items():
            distnum = int(dist[1:])
            for shape_name in shape_names:
                # this gets the shape predictions for all sequences with the
                # distance of 'dist'
                shapes = getattr(shape_dist, shape_name.lower())
                for seqid in shapes:
                    # Roll and HelT are calculated with two central bp steps
                    seqint = int(seqid)
                    cur = shapes[seqid]
                    b1 = self.bsites[0][seqint]
                    b2 = self.bsites[1][seqint] + 1
                    if seqint not in features:
                        features[seqint] = {}
                    #print(list(cur[b1-span_out:b1+span_in]),list(cur[b2-span_in:b2+span_out]))
                    b1_left = cur[b1-span_out:b1]
                    for i in range(0,len(b1_left)):
                        type = "%s_left_%d" % (shape_name,(i-span_out))
                        features[seqint][type] = b1_left[i]
                    b1_right = cur[b1:b1+span_in]
                    for i in range(0,span_in):
                        type = "%s_left_%d" % (shape_name,i)
                        features[seqint][type] = b1_right[i]

                    b2_left = cur[b2-span_in:b2]
                    for i in range(0,len(b2_left)):
                        type = "%s_right_%d" % (shape_name,(i-span_in))
                        features[seqint][type] = b2_left[i]
                    b2_right = cur[b2:b2+span_out]
                    for i in range(0,span_out):
                        type = "%s_right_%d" % (shape_name,i)
                        features[seqint][type] = b2_right[i]

        df_ret = pd.DataFrame.from_dict(features,orient='index')
        if as_tbl:
            return df_ret
        else:
            return df_ret.to_dict('records')
"""
