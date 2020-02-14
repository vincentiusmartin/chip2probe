import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import math
from matplotlib.backends.backend_pdf import PdfPages

import utils

class ProbeData:
    def __init__(self, probe, negctrl, sep="\t", percentile=95, toexp=False, idcol = None):
        """
        TODO: decide whether the file is in exp or log form
        IDCOL
        """
        self.cutoff = self.cutoff_from_negctrl(negctrl,percentile=percentile,sep=sep,toexp=toexp)

        if not isinstance(probe, pd.DataFrame):
            df = pd.read_csv(probe,delimiter=sep)
        else:
            df = probe.copy()

        cols2exp = ['Median', 'Median_o1', 'Median_o2', 'o1_r1', 'o1_r2', 'o1_r3',
           'o2_r1', 'o2_r2', 'o2_r3']
        if toexp:
            df = pd.concat([df[['Name']],np.exp(df[cols2exp]),df[['Sequence']]],axis=1)
        else:
            df = pd.concat([df[['Name']],df[cols2exp],df[['Sequence']]],axis=1)

        m1 = df[df['Name'].str.endswith("m1")].reset_index(drop=True)
        m1['Name'] = m1['Name'].apply(lambda x : x.rsplit("_",1)[0])
        m1.sort_values(by=['Name'], inplace=True)
        m1 = m1.reset_index(drop=True)

        m2 = df[df['Name'].str.endswith("m2")].reset_index(drop=True)
        m2['Name'] = m2['Name'].apply(lambda x : x.rsplit("_",1)[0])
        m2.sort_values(by=['Name'], inplace=True)
        m2 = m2.reset_index(drop=True)

        m3 = df[df['Name'].str.endswith("m3")].reset_index(drop=True)
        m3['Name'] = m3['Name'].apply(lambda x : x.rsplit("_",1)[0])
        m3.sort_values(by=['Name'], inplace=True)
        m3 = m3.reset_index(drop=True)

        wt = df[df['Name'].str.endswith("wt")].reset_index(drop=True)
        wt['Name'] = wt['Name'].apply(lambda x : x.rsplit("_",1)[0])
        wt.sort_values(by=['Name'], inplace=True)
        wt = wt.reset_index(drop=True)

        # Check if the order of every name column is the same
        zipped_names = zip(m1['Name'], m2['Name'], m3['Name'], wt['Name'])
        is_equal = all(len(set(x)) == 1 for x in zipped_names)
        if not is_equal:
            raise Exception("Name columns are different between wt, m1, m2, m3")

        self.table = {"m1":m1,"m2":m2,"m3":m3,"wt":wt}
        self.indivsum,self.twosites = self.make_replicas_permutation()
        self.medians = self.get_median_orientation()


    # ======== CONSTRUCTOR FUNCTION ========
    def cutoff_from_negctrl(self,negctrl,sep="\t",percentile=95,toexp=False):
        """
        get cutoff from negative control file
        negctrl: can be path or data frame
        """
        negcutoff = {}
        for orientation in [1,2]:
            colnames = ["o%d_r1" % orientation,"o%d_r2" % orientation,"o%d_r3" % orientation]
            if not isinstance(negctrl, pd.DataFrame):
                negdf = pd.read_csv(negctrl,delimiter=sep)
            else:
                negdf = negctrl.copy()
            negdf = negdf[colnames].values.tolist()
            if toexp:
                flatten = [np.exp(item) for sublist in negdf for item in sublist]
            else:
                flatten = [item for sublist in negdf for item in sublist]
            negcutoff["o%d"%orientation] = np.percentile(flatten,percentile)
            #h = sorted(flatten)
            #fit = stats.norm.pdf(h, np.mean(h), np.std(h))  #this is a fitting indeed
            #plt.plot(h,fit)
            #plt.axvline(negcutoff["o%d"%orientation],color='red',linestyle='dashed')
            #plt.show()
        return negcutoff

    def get_median_orientation(self):
        """
        Get only the median column from the table. We have median for orientation 1
        (o1) and orientation 2 (o2). Median from both which is treated as a mean
        value of o1 and o2 denoted as o0.

        return:
            dictionary of [xx]o[y] with xx=m1/m2/m3/wt and y=0/1/2
        """
        med_dict = {}
        orientations = [0,1,2]

        for orientation in orientations:
            if orientation == 0:
                colname = "Median"
            else:
                colname = "Median_o%d" % orientation
            med_dict["m1o%d"%orientation] = self.table['m1'][colname]
            med_dict["m2o%d"%orientation] = self.table['m2'][colname]
            med_dict["m3o%d"%orientation] = self.table['m3'][colname]
            med_dict["wto%d"%orientation] = self.table['wt'][colname]

        return med_dict

    def make_replicas_permutation(self):
        """
        Make permutation for each replicates, useful for hypothesis testing.
        Since individual sites is m1-m3+m2-m3 and two sites is wt - m3, we can take
        m3 from both and have: individual sites = m1 + m2 - m3 and two sites = wt.

        return:
            indivsum = permutations from m1,m2,m3
            twosites = permutations from wt
        """
        indivsum = {}
        twosites = {}

        for orientation in [1,2]:
            replica_cols = ["o%d_r1"%orientation,"o%d_r2"%orientation,"o%d_r3"%orientation]

            indivsum_permut = list(itertools.product(*[replica_cols]*3))
            twosites_permut = list(itertools.product(*[replica_cols]*1))

            indivsum_list = []
            for permut in indivsum_permut:
                indiv_sum = self.table["m1"][permut[0]] + self.table["m2"][permut[1]] - self.table["m3"][permut[2]]
                indivsum_list.append(indiv_sum.tolist())
            indivsum["o%d"%orientation] = pd.DataFrame(indivsum_list).transpose()

            twosites_list = []
            for permut in twosites_permut:
                twosites_sum = self.table["wt"][permut[0]]
                twosites_list.append(twosites_sum.tolist())
            twosites["o%d"%orientation] = pd.DataFrame(twosites_list).transpose()

        return indivsum,twosites

    # ======== ANALYSIS FUNCTION ========

    def check_cutoff_negctrl(self):
        print("Total sequences %d" % self.table["wt"].shape[0])
        for orientation in [1,2]:
            ori = "o%d"%orientation
            print("Orientation %s:" % ori)
            c = self.cutoff[ori]
            m1_med = self.medians["m1o%d"%orientation].tolist()
            m2_med = self.medians["m2o%d"%orientation].tolist()
            m3_med = self.medians["m3o%d"%orientation].tolist()
            wt_med = self.medians["wto%d"%orientation].tolist()
            m1_ct, m2_ct, wt_ct, m3_ct, total_ct = 0,0,0,0,0
            for i in range(len(m1_med)):
                if (m1_med[i] < c or m2_med[i] < c or wt_med[i] < c or m3_med[i] > c):
                    #print("  Sequence %d m1 %d m2 %d wt %d m3 %d" % (i, m1_med[i] < c, m2_med[i] < c, wt_med[i] < c, m3_med[i] > c))
                    if m1_med[i] < c:
                        m1_ct += 1
                    if m2_med[i] < c:
                        m2_ct += 1
                    if wt_med[i] < c:
                        wt_ct += 1
                    if m3_med[i] > c:
                        m3_ct += 1
                    total_ct += 1
            print("  Number of sequences violating cutoff m1 %d m2 %d wt %d m3 %d (total: %d)" % (m1_ct, m2_ct, wt_ct, m3_ct, total_ct))

    def get_seq(self,seqtype,indexes=[], othercols=False, tofile=False):
        """
        Get sequence of type 'seqtype' and a list containing the desired indexes.
        By default indexes is an emtpy list.
        Indexes can be set to -1 to return all sequences.
        othercols: information from other column is needed can be shown as well,
            if this is the case, the result will be a dictionary
        """
        if indexes == -1:
            indexes = self.table[seqtype].index.values
        seqlist = self.table[seqtype]['Sequence'][indexes].tolist()
        if not othercols:
            seqdict = {int(indexes[i]):seqlist[i] for i in range(len(seqlist))}
        else:
            other = self.table[seqtype][othercols].to_dict("list")
            seqdict = {}
            for i in range(len(seqlist)):
                content = {key:other[key][int(indexes[i])] for key in other}
                content["sequence"] = seqlist[i]
                seqdict[int(indexes[i])] = content
        #
        if not tofile:
            return seqdict
        else:
            if other:
                raise Exception("tofile feature is not compatible with othercols")
            keys = sorted(seqdict.keys())
            with open("sequences.txt",'w') as f:
                for key in keys:
                    f.write(">%s\n"%key)
                    f.write("%s\n"%seqdict[key])

    def get_mutpos(self,indexes=[]):
        wt = self.table['wt']['Sequence']
        m1 = self.table['m1']['Sequence']
        m2 = self.table['m2']['Sequence']

        if not indexes:
            indexes = wt.index.values

        mutposdict = {}
        for index in indexes:
            wtseq = wt[[index]].iloc[0]
            m1seq = m1[[index]].iloc[0]
            m2seq = m2[[index]].iloc[0]
            diff1 = [u+1 for u in range(len(wtseq)) if wtseq[u] != m1seq[u]]
            diff2 = [u+1 for u in range(len(wtseq)) if wtseq[u] != m2seq[u]]
            mutposdict[index] = diff1+diff2

        return mutposdict

    def scatter_boxplot_permutation(self,rownum):
        """
        make boxplot with scatterplot for individual vs twosites for a row
        """
        for orientation in [1,2]:
            # use loc of iloc as we want to access by index
            indivsum_df = self.indivsum["o%d"%orientation].loc[[rownum]].values[0]
            twosites_df = self.twosites["o%d"%orientation].loc[[rownum]].values[0]

            alldf = [indivsum_df,twosites_df]
            bp = plt.boxplot(alldf,positions = [1,1.5], widths=0.4)
            plt.xticks([1, 1.5], ['individual sum', 'two sites'])
            plt.setp(bp['boxes'], color='black')
            plt.setp(bp['caps'], color='black')

            for i in range(len(alldf)):
                y = alldf[i]
                x = np.random.normal(1+i*0.5, 0.02, size=len(y))
                plt.plot(x, y, 'r.', alpha=0.4,c='red')

            plotfilename = "row%so%d-box.png" % (rownum,orientation)
            print("Save distribution of row %s to %s" % (rownum,plotfilename))
            plt.savefig(plotfilename,positions=[0, 1])
            plt.clf() # clear canvas
        return plotfilename

    def scatter_boxplot(self, dflist, plotfilename = "box.png", log=False, ax=False):
        if log:
            alldf = [np.log(df) for df in dflist]
        else:
            alldf = dflist
        pos = np.linspace(1,1+len(alldf)*0.5-0.5,len(alldf))
        if not ax:
            bp = plt.boxplot(alldf, positions=pos, widths=0.4)
            plt.xticks(pos, ['wt_o1', 'm1_o1','m2_o1','m3_o1', 'wt_o2','m1_o2','m2_o2','m3_o2'])
            plt.setp(bp['boxes'], color='black')
            plt.setp(bp['caps'], color='black')
        else:
            bp = ax.boxplot(alldf, positions=pos, widths=0.4)
            ax.set_xticks(pos)
            ax.set_xticklabels(['wt_o1', 'm1_o1','m2_o1','m3_o1', 'wt_o2','m1_o2','m2_o2','m3_o2'])
            #ax.setp(bp['boxes'], color='black')
            #ax.setp(bp['caps'], color='black')

        for i in range(len(alldf)):
            y = alldf[i]
            x = np.random.normal(1+i*0.5, 0.02, size=len(y))
            if not ax:
                plt.plot(x, y, 'r.', alpha=0.4,c='red')
            else:
                ax.plot(x, y, 'r.', alpha=0.4,c='red')

        for i in range(len(alldf)):
            y = alldf[i]
            x = np.random.normal(1+i*0.5, 0.02, size=len(y))
            if not ax:
                plt.plot(x, y, 'r.', alpha=0.4,c='red')
            else:
                ax.plot(x, y, 'r.', alpha=0.4,c='red')

        if not ax:
            print("Save boxplot to %s" % plotfilename)
            plt.savefig(plotfilename,positions=[0, 1])
            plt.clf() # clear canvas
            return plotfilename
        else:
            # not necessary needed
            return ax

    def scatter_boxplot_row(self,rownum,log=False,ax=False):
        m1,m2,m3,wt = self.table['m1'], self.table['m2'], self.table['m3'], self.table['wt']
        m1_o1 = m1[['o1_r1','o1_r2','o1_r3']].loc[[rownum]].values[0]
        m2_o1 = m2[['o1_r1','o1_r2','o1_r3']].loc[[rownum]].values[0]
        m3_o1 = m3[['o1_r1','o1_r2','o1_r3']].loc[[rownum]].values[0]
        wt_o1 = wt[['o1_r1','o1_r2','o1_r3']].loc[[rownum]].values[0]
        m1_o2 = m1[['o2_r1','o2_r2','o2_r3']].loc[[rownum]].values[0]
        m2_o2 = m2[['o2_r1','o2_r2','o2_r3']].loc[[rownum]].values[0]
        m3_o2 = m3[['o2_r1','o2_r2','o2_r3']].loc[[rownum]].values[0]
        wt_o2 = wt[['o2_r1','o2_r2','o2_r3']].loc[[rownum]].values[0]

        #wt_seq = wt[['Sequence']].loc[rownum].values[0]
        #m1_seq = m1[['Sequence']].loc[rownum].values[0]
        #m2_seq = m2[['Sequence']].loc[rownum].values[0]
        #m3_seq = m3[['Sequence']].loc[rownum].values[0]

        #print("wt",wt_seq)
        #print("m1",m1_seq)
        #print("m2",m2_seq)
        #print("m3",m3_seq)

        alldf = [wt_o1,m1_o1,m2_o1,m3_o1, wt_o2,m1_o2,m2_o2,m3_o2]
        plotname = "row%s-box.png" % (rownum) # this is only used when ax=False
        return self.scatter_boxplot(alldf,plotfilename=plotname,log=log,ax=ax)

    def scatter_boxplot_distribution(self,rowlist,log=False,ax=False):
        m1,m2,m3,wt = self.table['m1'], self.table['m2'], self.table['m3'], self.table['wt']
        m1_o1 = m1[['o1_r1','o1_r2','o1_r3']].loc[rowlist].values.reshape(-1).tolist()
        m2_o1 = m2[['o1_r1','o1_r2','o1_r3']].loc[rowlist].values.reshape(-1).tolist()
        m3_o1 = m3[['o1_r1','o1_r2','o1_r3']].loc[rowlist].values.reshape(-1).tolist()
        wt_o1 = wt[['o1_r1','o1_r2','o1_r3']].loc[rowlist].values.reshape(-1).tolist()
        m1_o2 = m1[['o2_r1','o2_r2','o2_r3']].loc[rowlist].values.reshape(-1).tolist()
        m2_o2 = m2[['o2_r1','o2_r2','o2_r3']].loc[rowlist].values.reshape(-1).tolist()
        m3_o2 = m3[['o2_r1','o2_r2','o2_r3']].loc[rowlist].values.reshape(-1).tolist()
        wt_o2 = wt[['o2_r1','o2_r2','o2_r3']].loc[rowlist].values.reshape(-1).tolist()
        alldf = [wt_o1,m1_o1,m2_o1,m3_o1, wt_o2,m1_o2,m2_o2,m3_o2]
        plotname="bow-dist.png"
        return self.scatter_boxplot(alldf,plotfilename=plotname,log=log,ax=ax)


    def multi_scatter_boxplot(self, rowlist, log=False, filepath="scatter_boxplots.pdf"):
        """
        The first plot will always be the summary plot
        """
        numcol = 4
        numrow = 4
        fig, ax = plt.subplots(numrow, numcol, figsize=(25, 5))
        plt.subplots_adjust(hspace = 0.4, wspace=0.2)
        with PdfPages(filepath) as pdf:
            fig = plt.figure(figsize=(25,14))
            fig.subplots_adjust(hspace=0.4,wspace=0.2)
            n = 1
            cur_ax = fig.add_subplot(numcol,numrow,n)
            self.scatter_boxplot_distribution(rowlist, log=True, ax=cur_ax)
            for i in range(0,len(rowlist)):
                if n == 0:
                    fig = plt.figure(figsize=(25,14))
                    fig.subplots_adjust(hspace=0.4,wspace=0.2)
                n+=1
                cur_ax = fig.add_subplot(numcol,numrow,n)
                self.scatter_boxplot_row(rowlist[i],log=log,ax=cur_ax)
                cur_ax.set_title(rowlist[i])
                if n == numcol*numrow:
                    pdf.savefig(fig)
                    plt.close()
                    n = 0
            pdf.savefig(fig)
        plt.clf()
        plt.close()
