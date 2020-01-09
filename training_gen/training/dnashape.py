import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

from matplotlib.backends.backend_pdf import PdfPages

from teacup import utils
from teacup.training import trainingparser

class DNAShape:

    def __init__(self, path): #bpos1_list, bpos2_list
        # TODO: Checking in here
        for root, dirs, files in os.walk(path):
            for filename in files:
                path = "%s/%s" % (root,filename)
                filename, file_extension = os.path.splitext(path)
                if file_extension == ".ProT":
                    self.prot = utils.read_dictlist_file(path)
                elif file_extension == ".MGW":
                    self.mgw = utils.read_dictlist_file(path)
                elif file_extension == ".Roll":
                    self.roll = utils.read_dictlist_file(path)
                elif file_extension == ".HelT":
                    self.helt = utils.read_dictlist_file(path)

    def get_shape_names():
        return ["ProT","MGW","Roll","HelT"]

    def plot_average(self, labels, bsites, pthres=0.05, path=".", plotlabel="Average DNA shape"):
        plt.clf()
        colors = ['orangered','dodgerblue','lime']

        shapes = {"Propeller twist ":self.prot,"Helix twist":self.helt,"Roll":self.roll,"Minor Groove Width":self.mgw}
        #shapes = {"Propeller twist ":self.prot}
        keystolabel = {"additive":"*","cooperative":"*"} # ã ç

        min_dist = min(bsites[0].values())
        max_dist = max(bsites[0].values())
        # get the maximum seq length

        fig = plt.figure(figsize=(12,12))
        #with PdfPages("shape.pdf") as pdf:
        n = 0
        for sh in shapes: # shapes
            c = 0
            n += 1
            ax = fig.add_subplot(2,2,n)
            trimlen = len(next(iter(shapes[sh].values()))) - (max_dist - min_dist)
            yall = {}
            lshift_dict = {}
            curbsites = []

            # ==== Plotting part using 25, 50, and 75 percentiles ====
            # labels = cooperative or additive
            keys = []
            pr = True
            for label in labels:
                indexes = labels[label]
                curlist = []
                for key in indexes:
                    # we need to align the binding sites
                    # left shift is from the diff between first bsite and min bsite
                    lshift = bsites[0][key] - min_dist
                    lshift_dict[key] = lshift
                    aligned = shapes[sh][str(key)][lshift:]
                    aligned = aligned[:trimlen]
                    curlist.append(aligned)
                    if pr: # just to save the bpos location
                        curbsites = [min_dist, bsites[1][key] - lshift]
                        pr = False
                keys.append(label) # for marking significance
                yall[label] = curlist

                seqlen = len(curlist[0])
                if not all(len(l) == seqlen for l in curlist):
                    raise ValueError('not all lists have same length!')

                ylist = np.median(curlist, axis=0)
                xlist = [i+1 for i in range(0,seqlen)]
                y25p = np.percentile(curlist, 25, axis=0)
                y75p = np.percentile(curlist, 75, axis=0)
                ax.plot(xlist, ylist, alpha=0.8, label=label, c=colors[c], marker='o')
                ax.fill_between(xlist, y75p, y25p, alpha=0.15, facecolor=colors[c])
                c += 1

            # ==== Hypothesis testing to mark significant binding sites ====
            signiflabel = []
            for i in range(0,trimlen):
                # for now assume yall is of size 2
                arr_coop = [seq[i] for seq in yall['cooperative']]
                arr_add = [seq[i] for seq in yall['additive']]
                if not any(np.isnan(x) for x in arr_coop) and not any(np.isnan(x) for x in arr_add):
                    p1 = utils.wilcox_test(arr_coop,arr_add,"greater")
                    p2 = utils.wilcox_test(arr_coop,arr_add,"less")
                    if p1 <= pthres:
                        signiflabel.append(keystolabel["cooperative"])
                    elif p2 <= pthres:
                        signiflabel.append(keystolabel["additive"])
                    else:
                        signiflabel.append('')
                else:
                    signiflabel.append('')

            # ==== Mark binding sites as given from the input ====
            for m in curbsites:
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

# TODO FIX BASED ON THE CLASS

def plot_average_all(trainingpath,shapepath,distances):
    train = trainingparser.TrainingParser(trainingpath,motiflen=6)
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
