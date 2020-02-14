import rpy2.robjects as robjects
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from termcolor import colored
import statsmodels.stats.multitest as sm

#  custom library
import util.stats as sta

# ================ R function wrapper ================

def shapiro_test(df):
    pval_lists = df.apply(lambda row: sta.shapiro(sta.numeric(row.tolist())).rx("p.value")[0][0],axis=1)
    return pval_lists

def wilcox_test(df1,df2,wilcox_alternative='two.sided'):
    if df1.shape[0] != df2.shape[0]:
        print("Different number of rows!")
        return -1
    pval_lists = []
    for i in range(0,df1.shape[0]):
        row1 = sta.numeric(df1.loc[i].tolist())
        row2 = sta.numeric(df2.loc[i].tolist())
        pval_lists.append(sta.wilcox(row1,row2,alternative=wilcox_alternative)) # change to sta
    return pval_lists

# ================================

def wilcox_pthreshold(df1,df2,pval=0.05,alternative='two.sided'):
    wilcox = wilcox_test(df1,df2,alternative)
    p_series = pd.Series([True if p < pval else False for p in wilcox])
    return p_series

# ================ Internal for this program ================

def overlap(o1,o2):
    return sorted(list(set(o1)&set(o2)))

def const_satisfied(m1,m2,m3,wt,cutoff):
    constraint = (m1 > cutoff) &  (m2 > cutoff) & (wt > cutoff) & (m3 <= cutoff)
    return constraint

# ================ Tools ================

def classify_per_orientation(probes,pvalthres,classify_inconsistency=True, fdrcor=False):
    medians = probes.medians
    cutoff = probes.cutoff

    indivsum = probes.indivsum
    twosites = probes.twosites

    #pd.set_option('display.max_columns', None)
    #print(probes.table["wt"])
    #names = probes.table["wt"]["Name"].tolist() --> debugging purpose

    classification = {}
    pval = {}
    for orientation in [1,2]:
        ori = "o%d"%orientation
        m1_med = medians["m1o%d"%orientation].tolist()
        m2_med = medians["m2o%d"%orientation].tolist()
        m3_med = medians["m3o%d"%orientation].tolist()
        wt_med = medians["wto%d"%orientation].tolist()
        wilcox_p = wilcox_test(twosites[ori],indivsum[ori],'greater')
        if fdrcor:
            # index 1 is the corrected pval
            wilcox_p = sm.fdrcorrection(wilcox_p)[1]
        # previously done: if wilcox_p[i] < pvalthres and m1_med[i] > cutoff[ori] and m2_med[i] > cutoff[ori] and wt_med[i] > cutoff[ori] and m3_med[i] <= cutoff[ori] --> remove m1 and m2 for ets
        classification["cooperative_o%d"%orientation] = [i for i in range(0,len(wilcox_p)) if wilcox_p[i] < pvalthres and wt_med[i] > cutoff[ori] and m3_med[i] <= cutoff[ori]]
        print("Wilcox greater test orientation %d, # cooperative rows with p-val less than %.3f: %d/%d" \
                % (orientation,pvalthres,len(classification["cooperative_o%d"%orientation]),len(wilcox_p)))
        #print("COOPERATIVE\n",[names[i] for i in classification["cooperative_o%d"%orientation] if "dist" in names[i]])

        wilcox_p_less = wilcox_test(twosites[ori],indivsum[ori],'less')
        if fdrcor:
            wilcox_p_less =sm.fdrcorrection(wilcox_p_less)[1]
        classification["anticoop_o%d"%orientation] = [i for i in range(0,len(wilcox_p)) if wilcox_p_less[i] < pvalthres and wt_med[i] > cutoff[ori] and m3_med[i] <= cutoff[ori]]
        print("Wilcox less test orientation %d, # anticoop rows with p-val less than %.3f: %d/%d" \
                % (orientation,pvalthres,len(classification["anticoop_o%d"%orientation]),len(wilcox_p_less)))
        #print("\nANTICOOP\n",[names[i] for i in classification["anticoop_o%d"%orientation] if "dist" in names[i]])

        # get the rest as additive
        classification["additive_o%d"%orientation] = [i for i in range(0,indivsum[ori].shape[0]) if i not in classification["anticoop_o%d"%orientation] and i not in classification["cooperative_o%d"%orientation] and wt_med[i] > cutoff[ori] and m3_med[i] <= cutoff[ori]]
        print("Wilcox additive test orientation %d, # additive rows with p-val less than %.3f: %d/%d" \
                % (orientation,pvalthres,len(classification["additive_o%d"%orientation]),len(wilcox_p_less)))
        #print("\nADDITIVE\n",[names[i] for i in classification["additive_o%d"%orientation] if "dist" in names[i]])

        #leftover = [i for i in range(0,len(wilcox_p)) if i not in classification["cooperative_o%d"%orientation] and i not in classification["anticoop_o%d"%orientation] and i not in classification["additive_o%d"%orientation]]

    classification["cooperative_overlap"] = overlap(classification["cooperative_o1"],classification["cooperative_o2"])
    classification["anticoop_overlap"] = overlap(classification["anticoop_o1"],classification["anticoop_o2"])
    classification["additive_overlap"] = overlap(classification["additive_o1"],classification["additive_o2"])
    print("Number of overlap cooperative: %d" % len(classification["cooperative_overlap"]))
    print("Number of overlap anticoop: %d" % len(classification["anticoop_overlap"]))
    print("Number of overlap additive: %d" % len(classification["additive_overlap"]))

    if classify_inconsistency:
        classification["cooperative_o1_additive_o2"] = overlap(classification["cooperative_o1"],classification["additive_o2"])
        classification["cooperative_o2_additive_o1"] = overlap(classification["cooperative_o2"],classification["additive_o1"])

        classification["additive_o1_anticoop_o2"] = overlap(classification["additive_o1"],classification["anticoop_o2"])
        classification["additive_o2_anticoop_o1"] = overlap(classification["additive_o2"],classification["anticoop_o1"])

        classification["cooperative_o1_anticoop_o2"] = overlap(classification["cooperative_o1"],classification["anticoop_o2"])
        classification["cooperative_o2_anticoop_o1"] = overlap(classification["cooperative_o2"],classification["anticoop_o1"])

        print("Number of ooperative_o1_additive_o2: %d" % len(classification["cooperative_o1_additive_o2"]))
        print("Number of cooperative_o2_additive_o1: %d" % len(classification["cooperative_o2_additive_o1"]))
        print("Number of additive_o1_anticoop_o2: %d" % len(classification["additive_o1_anticoop_o2"]))
        print("Number of additive_o2_anticoop_o1: %d" % len(classification["additive_o2_anticoop_o1"]))
        print("Number of cooperative_o1_anticoop_o2: %d" % len(classification["cooperative_o1_anticoop_o2"]))
        print("Number of cooperative_o2_anticoop_o1: %d" % len(classification["cooperative_o2_anticoop_o1"]))

    return classification

# printlvl: WARN,ALL
# TODO: need edit
def classify_orientation_combined(probes,pvalthres,printlvl="WARN"):
    tbl = probes.table
    cutoff = probes.cutoff
    med = probes.medians
    indivsum = probes.indivsum
    twosites = probes.twosites

    m1,m2,m3,wt = tbl["m1"],tbl["m2"],tbl["m3"],tbl["wt"]
    numrow = wt.shape[0]

    m1_med_o1,m2_med_o1,m3_med_o1,wt_med_o1 = med["m1o1"],med["m2o1"],med["m3o1"],med["wto1"]
    const_o1 = const_satisfied(m1_med_o1,m2_med_o1,m3_med_o1,wt_med_o1,cutoff["o1"])
    coop_o1 = wilcox_pthreshold(twosites["o1"],indivsum["o1"],pvalthres,'greater')
    coop_o1 = const_o1 & coop_o1 # must satisfy constraint
    steric_o1 = wilcox_pthreshold(twosites["o1"],indivsum["o1"],pvalthres,'less')
    steric_o1 = const_o1 & steric_o1
    additive_o1 = const_o1 & ~coop_o1 & ~steric_o1

    m1_med_o2,m2_med_o2,m3_med_o2,wt_med_o2 = med["m1o2"],med["m2o2"],med["m3o2"],med["wto2"]
    const_o2 = const_satisfied(m1_med_o2,m2_med_o2,m3_med_o2,wt_med_o2,cutoff["o2"])
    coop_o2 =  wilcox_pthreshold(twosites["o2"],indivsum["o2"],pvalthres,'greater')
    coop_o2 = const_o2 & coop_o2 # must satisfy constraint
    steric_o2 =  wilcox_pthreshold(twosites["o2"],indivsum["o2"],pvalthres,'less')
    steric_o2 = const_o2 & steric_o2
    additive_o2 = const_o2 & ~coop_o2 & ~steric_o2

    print("Pass cutoff in both: %d" % sum(const_o1[i] and const_o2[i] for i in range(1,numrow+1)))

    warns = {"coopsteric":[],"notconstraint":[]}
    # CHOOSE MODEL TO USE
    training = [["rownum","indivsites","twosites","label"]]

    for i in range(1,numrow+1):
        twosites_o1 = wt_med_o1[i] - m3_med_o1[i]
        indivsum_o1 = m1_med_o1[i] + m2_med_o1[i] - 2*m3_med_o1[i]
        twosites_o2 =  wt_med_o2[i] - m3_med_o2[i]
        indivsum_o2 = m1_med_o2[i] + m2_med_o2[i] - 2*m3_med_o2[i]
        # Cooperative:
        if coop_o1[i] and coop_o2[i]:
            # we want to incorporate bigger range, so take the max
            chosenmax = max(twosites_o1,twosites_o2)
            if twosites_o1 == chosenmax:
                training.append([i+1,indivsum_o1,twosites_o1,"cooperative"])
            else: # twosites_o2 == chosenmax:
                training.append([i+1,indivsum_o2,twosites_o2,"cooperative"])
        elif coop_o1[i]:
            if steric_o2[i]:
                print("%s O1 is cooperative but o2 is steric for probe %d"%(colored('WARN:', 'red'),i))
                warns['coopsteric'].append(i)
            elif not const_o2[i] or additive_o2[i]:
                training.append([i,indivsum_o1,twosites_o1,"cooperative"])
                if printlvl=="ALL" and additive_o2:
                    print("O1 is cooperative but o2 is additive for probe %d"%i)
        elif coop_o2[i]:
            if steric_o1[i]:
                print("%s O2 is cooperative but o1 is steric for probe %d"%(colored('WARN:', 'red'),i))
                warns['coopsteric'].append(i)
            elif not const_o1[i] or additive_o1[i]:
                training.append([i,indivsum_o2,twosites_o2,"cooperative"])
                if printlvl=="ALL" and additive_o1:
                    print("O2 is cooperative but o1 is additive for probe %d"%i)
        # =============STERIC=============
        elif steric_o1[i] and steric_o2[i]:
            chosenmin = min(twosites_o1,twosites_o2)
            if twosites_o1 == chosenmin:
                training.append([i,indivsum_o1,twosites_o1,"steric"])
            else: # twosites_o2 == chosenmax:
                training.append([i,indivsum_o2,twosites_o2,"steric"])
        elif steric_o1[i]:
            if not const_o2[i] or additive_o2[i]:
                training.append([i,indivsum_o1,twosites_o1,"steric"])
                if printlvl=="ALL" and additive_o2[i]:
                    print("O1 is steric but o2 is additive for probe %d"%i)
            #if coop_o2[i]: -- not needed as it was covered on COOP part
            #    print("%s O1 is steric but o2 is cooperative for probe %d"%(colored('WARN:', 'red'),(i+1)))
        elif steric_o2[i]:
            if not const_o1[i] or additive_o1[i]:
                training.append([i+1,indivsum_o2,twosites_o2,"steric"])
                if printlvl=="ALL" and additive_o1[i]:
                    print("O2 is steric but o1 is additive for probe %d"%i)
            #if coop_o1[i]: -- not needed for the same reason as above
            #    print("%s O2 is steric but o1 is coperative for probe %d"%(colored('WARN:', 'red'),(i+1)))
        # =============ADDITIVE=============
        elif additive_o1[i] and additive_o2[i]:
            chosenmax = max(twosites_o1,twosites_o2)
            if twosites_o1 == chosenmax:
                training.append([i,indivsum_o1,twosites_o1,"additive"])
            else: # twosites_o2 == chosenmax:
                training.append([i,indivsum_o2,twosites_o2,"additive"])
        elif additive_o1[i] and not const_o2[i]:
            training.append([i,indivsum_o1,twosites_o1,"additive"])
        elif additive_o2[i] and not const_o1[i]:
            training.append([i,indivsum_o2,twosites_o2,"additive"])
        # ==============OTHERS==============
        elif not const_o1[i] and not const_o2[i]:
            warns['notconstraint'].append(i)
            if printlvl=="ALL":
                print("O1 and O2 don't meet constraint for probe %d, don't include"%(i))
        else:
            print("%s: Something's weird in row %d"%(colored('WARN:', 'red'),i))
            print("  O1: pass constraint: %s,cooperative: %s, steric: %s, additive: %s"%(const_o1[i],coop_o1[i],steric_o1[i],additive_o1[i]))
            print("  O2: pass constraint: %s,cooperative: %s, steric: %s, additive: %s"%(const_o2[i],coop_o2[i],steric_o2[i],additive_o2[i]))
            break
    print("Cooperative but steric count: %d" % len(warns['coopsteric']))
    print("Rows that don't satisfy constraint %d" % len(warns['notconstraint']))

    classification = {}
    classification["cooperative"] = [str(x[0]) for x in training if x[3] == "cooperative"]
    classification["additive"] = [str(x[0]) for x in training if x[3] == "additive"]
    classification["steric"] = [str(x[0]) for x in training if x[3] == "steric"]
    classification["coop_steric"] = [str(x) for x in warns['coopsteric']]

    return classification

# ================ PLOTTING ================

def plot_median_binding_sum(probes,classification,orientation,
                            log=True,
                            plotname="plot",
                            tfname="",
                            plotnonsignif=True,
                            mark=[],
                            subset_additive=[],
                            subset_coop=[],
                            subset_steric=[]): # dictionary with subset of additive/coop/steric
    """
    orientation: 0 is overlap, else the value (1 or 2 is used)
    """
    med = probes.medians
    m1,m2,m3,wt = med["m1o%d"%orientation],med["m2o%d"%orientation],med["m3o%d"%orientation],med["wto%d"%orientation]
    ori = "overlap" if orientation == 0 else "o%d"%orientation

    if subset_coop:
        coop_idxs = [int(x) for x in subset_coop]
    else:
        coop_idxs = [int(x) for x in classification["cooperative_%s"%ori]]

    if subset_additive:
        steric_idxs = [int(x)for x in subset_steric]
    else:
        steric_idxs = [int(x) for x in classification["anticoop_%s"%ori]]

    if subset_steric:
        additive_idxs = [int(x) for x in subset_additive]
    else:
        additive_idxs = [int(x) for x in classification["additive_%s"%ori]]

    # negative val?
    x_axis = (m1+m2-2*m3).tolist()
    y_axis = (wt-m3).tolist()

    if log:
        x_axis = np.log(x_axis)
        y_axis = np.log(y_axis)

    for i in coop_idxs:
        if x_axis[i] > y_axis[i]:
            print("-----")
            print("index",i,"additive",x_axis[i], "cooperative", y_axis[i])
            print(m1.loc[i], m2.loc[i], m3.loc[i])
            print("indiv",probes.indivsum[ori].loc[i].tolist(),"two sites",probes.twosites[ori].loc[i].tolist())
    x_coop = np.take(x_axis,coop_idxs)
    y_coop = np.take(y_axis,coop_idxs)
    x_steric = np.take(x_axis,steric_idxs)
    y_steric = np.take(y_axis,steric_idxs)
    x_additive = np.take(x_axis,additive_idxs)
    y_additive = np.take(y_axis,additive_idxs)

    mark_convert = [x-1 for x in mark]
    x_mark = np.take(x_axis,mark_convert)
    y_mark = np.take(y_axis,mark_convert)

    if plotnonsignif:
        todelete_idxs = steric_idxs+coop_idxs+additive_idxs
        x_leftover = np.delete(x_axis,todelete_idxs)
        y_leftover = np.delete(y_axis,todelete_idxs)
        plt.scatter(x_leftover,y_leftover,s=3,c='yellow',zorder=1) # s = 0.15 # yellow

    """
    Wilcox greater test orientation 1, # cooperative rows with p-val less than 0.050: 918/6708
    Wilcox less test orientation 1, # steric rows with p-val less than 0.050: 376/6708
    Wilcox greater test orientation 2, # cooperative rows with p-val less than 0.050: 919/6708
    Wilcox less test orientation 2, # steric rows with p-val less than 0.050: 397/6708
    Number of overlap coop: 137
    Number of overlap steric: 135
    Number of overlap additive: 575
    """

    # === Scatter ====
    plt.scatter(x_coop,y_coop,s=2,c='blue',label='coooperative') # xkcd:azure
    plt.scatter(x_additive,y_additive,s=2,c='gray',label='additive') # gray
    plt.scatter(x_steric,y_steric,s=2,c='red',label='anti-cooperative') # red
    plt.scatter(x_mark,y_mark,s=10,c='orange') # orange

    #print("#coop %d, #additive %d, #steric %d, #filtered %d, #total %d" % (len(x_coop),len(x_additive),len(x_steric),len(x_leftover),sum([len(x_coop),len(x_additive),len(x_steric),len(x_leftover)])))

    # Custom range:
    #plt.xlim(6,12)
    #plt.ylim(6,12)

    lims = [
        np.min([plt.xlim(), plt.ylim()]),  # min of both axes
        np.max([plt.xlim(), plt.ylim()]),  # max of both axes
    ]
    #plt.xlim((5,12))
    #plt.ylim((5,12))
    plt.plot(lims, lims, 'k-', alpha=0.75, c='black', linewidth=0.5) #zorder=0
    plt.legend(loc="upper left",markerscale=3)
    #---------------
    # now plot both limits against each other
    plt.xlabel('log(m1-m3+m2-m3)' if log else 'm1-m3+m2-m3', fontsize=13)
    plt.ylabel('log(wt-m3)' if log else 'wt-m3',fontsize=13)
    plt.title("Cluster binding for TF %s" % tfname)
    print("Save %s-scatter.png to file" % (plotname))
    plt.savefig(plotname+"-scatter.png",bbox_inches='tight')
    plt.clf() # clear canvas

    return pd.DataFrame({"individual sum":x_axis,"two sites":y_axis})
