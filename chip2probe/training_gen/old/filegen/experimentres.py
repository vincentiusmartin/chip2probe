'''
Created on Oct 9, 2019

@author: vincentiusmartin
'''

import pandas as pd
import numpy as np

import sys
sys.path.append("..")
# for bio
import utils

import matplotlib.pyplot as plt

def plot_chamber_corr(filepaths, id_prefix = "Coop1Ets", log=False):
    dfs = [pd.read_csv(path, sep="\t") for path in filepaths]
    dfs = [df[df["ID"].str.startswith(id_prefix).fillna(False)][['ID','Alexa488Adjusted']].reset_index(drop=True) for df in dfs]

    combs = [(i1, i2) for i1 in range(0,len(dfs) - 1) for i2 in range(i1 + 1,len(dfs)) if i1 != i2]
    for c in combs:
        i1 = c[0]
        i2 = c[1]
        df_cham = dfs[i1].set_index('ID').join(dfs[i2].set_index('ID'),lsuffix='_1', rsuffix='_2')
        if log:
            df_cham['Alexa488Adjusted_1'] = np.log(df_cham['Alexa488Adjusted_1'])
            df_cham['Alexa488Adjusted_2'] = np.log(df_cham['Alexa488Adjusted_2'])
        plt.scatter(df_cham['Alexa488Adjusted_1'], df_cham['Alexa488Adjusted_2'],s=1, c='blue')
        corr = df_cham['Alexa488Adjusted_1'].corr(df_cham['Alexa488Adjusted_2'])

        #plt.xlabel('Chamber-%d' % (i1 + 1))
        #plt.ylabel('Chamber-%d' % (i2 + 1))
        plt.xlabel('FL 150nM scan 550_5')
        plt.ylabel('FL 30nM scan 550_50')
        plt.title("Experiment correlattion")
        lims = [
            np.min([plt.xlim(), plt.ylim()]),  # min of both axes
            np.max([plt.xlim(), plt.ylim()]),  # max of both axes
        ]
        plt.plot(lims, lims, 'k-', alpha=2, c='red', linewidth=0.5) #zorder=0
        plt.text(3, plt.ylim()[1] - 0.1 * plt.ylim()[1], s = "R² = %f" % (corr*corr), fontsize=12)
        figname = "corr_log-%d-%d.png" % (i1+1,i2+1) if log else "corr-%d-%d.png" % (i1+1,i2+1)
        plt.savefig(figname)
        plt.clf()

def count_sim_strs(a,b):
    u=zip(a,b)
    x = 0
    for i,j in u:
        if i==j:
            x += 1
    return x


# note: ets1_A549_seq26_dist_5
def make_coopfile(filepath, id_prefix = "Coop1Ets", seqlen = 36, numrep = 3, log=False):
    df = pd.read_csv(filepath, sep="\t")
    ###

    dpref = df[df["ID"].str.startswith(id_prefix).fillna(False)][['ID','Name','Sequence','Alexa488Adjusted']]

    if log:
        dpref['Alexa488Adjusted'] = np.log(dpref['Alexa488Adjusted'])

    dpref["Group"] = dpref.apply (lambda row: "_".join(row["Name"].split("_")[:-3]), axis=1)
    g = dpref.groupby("Group")

    reslist = []
    nctrl_list = []
    i = 0
    for name, group  in g:
        keyid = {}
        gdict = {}
        """if len(group) != 24:
            print(name,group)
            print(len(group))
        """
        if i % 100 == 0:
            print("Progress: %d/%d %.2f%%" % (i,len(g),float(i)/len(g) * 100))
        i += 1
        for idx, row in group.iterrows():
            cur_name, mut_type, rep, ori = row["Name"].rsplit('_', 3)
            # Checking here for duplicated entry that's why we group per sequence
            curseq = row["Sequence"][:seqlen]
            revcurseq = utils.revcompstr(curseq)

            curkey = ""
            for key in gdict.keys():
                # check similarity with similarity score of 28
                if count_sim_strs(curseq,key) > 28 or count_sim_strs(revcurseq,key) > 28:
                    curkey = key
                    break
            if not curkey:
                curkey = curseq
                keyid[curkey] = len(gdict)
                gdict[curkey] = {}
            if mut_type not in gdict[curkey]:
                if keyid[curkey] > 0:
                    gdict[curkey][mut_type] = {"Name": "%s_dup%d_%s" % (cur_name,keyid[curkey], mut_type)}
                else:
                    gdict[curkey][mut_type] = {"Name": "%s_%s" % (cur_name, mut_type)}
            if "Sequence" not in gdict[curkey] and ori == "o1":
                # we save the positive strand here
                gdict[curkey][mut_type]["Sequence"] = row["Sequence"][:seqlen]
            gdict[curkey][mut_type]["%s_%s" % (ori,rep)] = row["Alexa488Adjusted"]
        for s in gdict:
            for gelm in gdict[s]:
                o1_intensity = [gdict[s][gelm]["o1_r%d" % i] for i in range(1,numrep+1)]
                o2_intensity = [gdict[s][gelm]["o2_r%d" % i]  for i in range(1,numrep+1)]
                gdict[s][gelm]["Median"] = np.median(o1_intensity + o2_intensity)
                gdict[s][gelm]["Median_o1"] = np.median(o1_intensity)
                gdict[s][gelm]["Median_o2"] = np.median(o2_intensity)

        for k in gdict.keys():
            for mtype in gdict[k]:
                name = gdict[k][mtype]["Name"]
                if "NegativeCtrl" in name:
                    nctrl_list.append(gdict[k][mtype])
                else:
                    reslist.append(gdict[k][mtype])
    print("Progress: %d/%d %.2f%%" % (i,len(g),float(i)/len(g) * 100))

    o_r_list =  ["%s_r%s" % (o,r) for o in ["o1","o2"] for r in range(1,numrep + 1)]
    res_df = pd.DataFrame(reslist, columns=["Name","Median","Median_o1","Median_o2"] + o_r_list + ["Sequence"]).sort_values(by=['Name'])
    nctrl_df = pd.DataFrame(nctrl_list, columns=["Name","Median","Median_o1","Median_o2"] + o_r_list + ["Sequence"]).sort_values(by=['Name'])
    return res_df,nctrl_df
