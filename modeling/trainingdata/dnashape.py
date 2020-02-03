import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

import util.util as util
import util.stats as st

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

class DNAShape:

    def __init__(self, path): #bpos1_list, bpos2_list
        # TODO: Checking in here
        for root, dirs, files in os.walk(path):
            for filename in files:
                path = "%s/%s" % (root,filename)
                filename, file_extension = os.path.splitext(path)
                if file_extension == ".ProT":
                    self.prot = util.read_dictlist_file(path)
                elif file_extension == ".MGW":
                    self.mgw = util.read_dictlist_file(path)
                elif file_extension == ".Roll":
                    self.roll = util.read_dictlist_file(path)
                elif file_extension == ".HelT":
                    self.helt = util.read_dictlist_file(path)

    def get_shape_names():
        return ["ProT","MGW","Roll","HelT"]

    def get_sequence_count(self, seqlist, avg=True):
        it = iter(seqlist)
        seqlen = len(next(it))
        if not all(len(l) == seqlen for l in it):
            raise ValueError('not all lists have same length!')
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
    def plot_average(self, labels, site1list, site2list, sequences, pthres=0.05, path="shape.pdf", plotlabel="Average DNA shape"):
        plt.clf()
        #colors = ['orangered','dodgerblue','lime']
        s1 = int(sum(site1list.values()) / len(site1list.values()))
        s2 = int(sum(site2list.values()) / len(site2list.values()))

        shapes = {"Minor Groove Width":self.mgw}
        keystolabel = {"additive":"ã","cooperative":"ç"} # ã ç

        if "cooperative" in labels:
            cooperative_list = [sequences[idx] for idx in labels["cooperative"]]
            nuc_ct_coop = self.get_sequence_count(cooperative_list, avg=False) #self.get_sequence_count(list(sequences.values()), avg=True)
        if "additive" in labels:
            additive_list = [sequences[idx] for idx in labels["additive"]]
            nuc_ct_add = self.get_sequence_count(additive_list, avg=False)


        # get the maximum seq length
        fig = plt.figure(figsize=(12,12))
        n = 0 # for the subplot
        colors = {"additive":'orangered',"cooperative":'dodgerblue'}
        for sh in shapes: # shapes
            n += 1
            ax = fig.add_subplot(2,2,n)
            yall = {}
            lshift_dict = {}
            flen = int(np.mean([len(x) for x in shapes[sh].values()]))

            # ==== Plotting part using 25, 50, and 75 percentiles ====
            # labels = cooperative or additive
            miny = float("inf")
            for label in labels:
                indexes = labels[label]
                curlist = []
                for key in indexes:
                    shape = shapes[sh][str(key+1)]
                    curlist.append(shape)
                yall[label] = curlist
                seqlen = len(curlist[0])
                if not all(len(l) == seqlen for l in curlist):
                    raise ValueError('not all lists have same length!')

                ylist = np.median(curlist, axis=0)
                xlist = [i+1 for i in range(0,seqlen)]
                y25p = np.percentile(curlist, 25, axis=0)
                y75p = np.percentile(curlist, 75, axis=0)
                ax.plot(xlist, ylist, alpha=0.8, label=label, c=colors[label], marker='o')
                ax.fill_between(xlist, y75p, y25p, alpha=0.15, facecolor=colors[label])

            # ==== Hypothesis testing to mark significant binding sites ====
            signiflabel = []
            if "cooperative" in labels and "additive" in labels:
                for i in range(0,flen):
                    # for now assume yall is of size 2
                    arr_coop = [seq[i] for seq in yall['cooperative']]
                    arr_add = [seq[i] for seq in yall['additive']]
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

            # ==== Mark binding sites as given from the input ====
            for m in [s1,s2]:
                ax.axvline(x=m,linewidth=1, color='g',linestyle='--')

            ax.set_xlim(1,seqlen)
            xi = [x for x in range(1,seqlen+1)]
            #label = ['' if x not in signiflist else '*' for x in xi]
            ax.set_xticks(xi)
            ax.set_xticklabels(signiflabel)
            ax.yaxis.set_label_text(sh)
            ax.xaxis.set_label_text('Sequence position')
            ax.legend(loc="upper right")
            ax.set_title(plotlabel)
            low_y = ax.get_ylim()[0]
            hfactor = ax.get_ylim()[1] -  ax.get_ylim()[0]
            ax.set_ylim(bottom = low_y - np.abs(0.35 * hfactor))
            for i in range(0,4): # low_y here?
                ax.yaxis.get_major_ticks()[i].label1.set_visible(False)

            #for ymaj in ax.yaxis.get_major_ticks():
            #    print(ymaj.label2)
            #    ymaj.label1.set_visible(False)

            base_color = {'A': "red", 'C': "green" , 'G': "orange", 'T': "blue"}
            bot_anchor = 0

            nuc_dict = {}
            if "cooperative" in labels:
                nuc_dict["co"] = nuc_ct_coop
            if "additive" in labels:
                nuc_dict["ad"] = nuc_ct_add

            #lv = []
            #for key in nuc_dict:
            #    for i in range(0,len(nuc_dict[key])):
            #        lv.extend(list(nuc_dict["co"][0].values()))
            maxct = len(site1list)

            for base in base_color:
                for key in nuc_dict:
                    nuc_ct = nuc_dict[key]
                    inset_ax = inset_axes(ax,
                                  height="5%",
                                  width="100%",
                                  bbox_to_anchor= (0, bot_anchor, 1, 1),
                                  bbox_transform=ax.transAxes,
                                  loc=8)
                    inl = [elm[base] for elm in nuc_ct]
                    # also add 0.5 here to center
                    xbar = [i+0.5 for i in range(0,len(inl))] # should use the same axis with ax but this works...
                    sns.barplot(x = xbar, y = inl, color = base_color[base], ax=inset_ax)
                    inset_ax.get_xaxis().set_visible(False)
                    inset_ax.set_ylim(top=maxct)
                    inset_ax.set_yticklabels([])
                    inset_ax.patch.set_visible(False)
                    inset_ax.set_ylabel("%s_%s"%(base,key),rotation=0)
                    inset_ax.yaxis.label.set_color(base_color[base])
                    bot_anchor += 0.05
            # (A, green), thymine (T, red), cytosine (C, orange), and guanine (G, blue).

        with PdfPages(path) as pdf:
            pdf.savefig(fig)

# ============== This is a separate class to contain everything ==============

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
