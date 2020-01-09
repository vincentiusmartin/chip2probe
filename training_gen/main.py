'''
Created on Oct 9, 2019

@author: vincentiusmartin
'''
import sys
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe") # PATH TO UTIL

from filegen.experimentres import plot_chamber_corr, make_coopfile
from probes import classifier
import utils
from probes.probedata import ProbeData
import random

import pandas as pd
import re, os

# TODO: make this more general
def filtertrain_csv(probes, categories, pfiles):
    seqdict =  probes.get_seq("wt",indexes=categories["coop_overlap"])

    for p in pfiles:
        d = pd.read_csv(p)
        print(d.columns)

# use this if we have any duplicated colunnbs
# Attr keep=False: Mark all duplicates as True, else the first occurence
# won't be considered.
#dup = feature_df[feature_df.duplicated(['sequence'], keep=False)]
#pd.set_option('display.max_columns', 500)
#dup.to_csv("seqdata.csv")
def seqdata_to_feature_dict(filepaths):
    dflist1 = []
    for fpath in filepaths:
        with open(fpath) as f:
            cols = f.readline()
        if "\t" in cols:
            sep = "\t"
        elif "," in cols:
            sep = ","
        df = pd.read_csv(fpath, sep=sep)
        all_col = {"wt", "core1_start", "core1_end", "site1_pref", "core2_start", "core2_end",  "site2_pref", "distance"}.issubset(set(df.columns)) # .issubset
        if all_col:
            df = df[["wt", "core1_start", "core1_end", "site1_pref", "core2_start", "core2_end",  "site2_pref", "distance"]]
            dflist1.append(df)
        else:
            print("not all columns found in %s" % fpath)

    feature_df = pd.concat(dflist1)
    #feature_df.to_csv("feature_test.tsv",sep="\t")
    feature_df = feature_df.rename(columns={"wt": "sequence"})
    nodup_feature_df = feature_df.groupby('sequence').median().reset_index()
    nodup_feature_df['site1_pos'] = nodup_feature_df.apply(lambda row: int((row["core1_start"] + row["core1_end"]) // 2), axis=1)
    nodup_feature_df['site2_pos'] = nodup_feature_df.apply(lambda row: int((row["core2_start"] + row["core2_end"]) // 2), axis=1)
    nodup_feature_df = nodup_feature_df[["sequence", "site1_pos", "site2_pos", "site1_pref", "site2_pref", "distance"]]

    # transform to dictionary
    fdict = nodup_feature_df.set_index('sequence').T.to_dict('dict')
    return fdict

def make_training(analysis_respath, seqwithin_pattern, probedata, classification, print_not_found=False, corelen = 4):
    # core len is hardcoded , TODO: fix later
    seqwithin_re_pattern = re.compile(seqwithin_pattern)
    # this first to get imads pref, bsite pos, and bsite distance--i.e. all training features
    seqwithin_paths = []
    for path, dirs, filenames in os.walk(analysis_respath):
        for fname in filter(lambda name: seqwithin_re_pattern.match(name),filenames):
            filepath = os.path.join(path, fname)
            seqwithin_paths.append(filepath)
    # if the directory contains the pattern we want
    if seqwithin_paths:
        fdict = seqdata_to_feature_dict(seqwithin_paths) # sep="\t"
        # TODO: put if not...

    # now we get sequence and all the features, we now label them using classification result
    training_data = []
    notfound_count = 0
    for key in classification:
        #if key.endswith("overlap"):#key.endswith("o1") or key.endswith("o2"):
        # we label cooperative_o1_anticoop_o2 as cooperative, so we can do this:
        label = key.split("_")[0]
        #pd.set_option('display.max_columns', 500)
        #print(probedata.table["wt"])

        seqs = probedata.get_seq("wt",classification[key], othercols=["Name"])
        for idx in seqs:
            if seqs[idx]["sequence"] in fdict:
                #if not key.endswith("o2"):
                curseq = seqs[idx]["sequence"]
                features = fdict[curseq]
                seqlen = len(curseq)
                # site 1 is the stronger site
                if fdict[curseq]['site1_pref'] > fdict[curseq]['site2_pref']:
                    features = {
                        'site_str_pos':fdict[curseq]['site1_pos'],
                        'site_wk_pos':fdict[curseq]['site2_pos'],
                        'site_str_score':fdict[curseq]['site1_pref'],
                        'site_wk_score':fdict[curseq]['site2_pref'],
                        'distance':fdict[curseq]['distance']
                    }
                else:
                    features = {
                        'site_str_pos':fdict[curseq]['site2_pos'],
                        'site_wk_pos':fdict[curseq]['site1_pos'],
                        'site_str_score':fdict[curseq]['site2_pref'],
                        'site_wk_score':fdict[curseq]['site1_pref'],
                        'distance':fdict[curseq]['distance']
                    }
                """
                ### we do overlap so this is not needed
                else: #if o2, reverse complement the sequence and flip the features
                    s = seqs[idx]["sequence"]
                    curseq = utils.revcompstr(s)
                    seqlen = len(curseq)
                    features = {
                        'site1_pos':seqlen-fdict[s]['site2_pos'],
                        'site2_pos':seqlen-fdict[s]['site1_pos'],
                        'site1_pref':fdict[s]['site2_pref'],
                        'site2_pref':fdict[s]['site1_pref'],
                        'distance':fdict[s]['distance']
                    }
                """
                ftrs = {**{'name':seqs[idx]["Name"],"sequence":curseq,"label":label, "id":"%s" % (idx)}, **features}
                training_data.append(ftrs)
            else: # if not found
                #if not "weak" in seqs[idx]["Name"] and not "dist"  in seqs[idx]["Name"]:
                #print("Not found",seqs[idx]["Name"])
                notfound_count += 1
                if print_not_found:
                    print("couldn't find %s in feature dict" % seqs[idx])

    training_df = pd.DataFrame(training_data, columns=["id","name","sequence", "site_str_pos", "site_wk_pos", "site_str_score", "site_wk_score", "distance", "label"])
    print("Total training rows: %d, not found count: %d" % (len(training_data),notfound_count))
    training_df.drop_duplicates(subset=['sequence'],keep='first', inplace=True)
    print("Number of row in training df after dropping duplicates: %d"%training_df.shape[0])

    training_df["site_str_pos"] = training_df["site_str_pos"].astype(int)
    training_df["site_wk_pos"] = training_df["site_wk_pos"].astype(int)
    training_df["distance"] = training_df["distance"].astype(int)
    #training_df.to_csv("training.tsv", sep="\t", index=False)
    return training_df

if __name__ == '__main__':
    dir = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191004_coop-PBM_Ets1_v1_1st/2.processed_gpr"
    filepaths = ["%s/%s" % (dir,"20191004_258614510001_ETS1_550_5_1-4_alldata.txt"), "%s/%s" % (dir,"20191004_258614510001_ETS1_550_5_2-4_alldata.txt"), "%s/%s" % (dir,"20191004_258614510001_ETS1_550_5_3-4_alldata.txt"), "%s/%s" % (dir,"20191004_258614510001_ETS1_660_60_4-4_alldata.txt")]
    probe_analysis_path = "/Users/vincentiusmartin/Research/chip2gcPBM/result"
    negcutoff = 95
    pvalthres = .01
    tfname = "Ets1"
    #filepaths = ["/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/processed_gpr/20191030_258614510001_550_50_ch1.gpr_alldata.txt"]
    #plot_chamber_corr(filepaths, log=False)

    infile = "/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/2.processed_gpr/20191004_258614510001_ETS1_550_5_1-4_alldata.txt"
    fname = os.path.basename(infile)
    fname_woext = os.path.splitext(fname)[0]

    print("Making multisites file for %s..." % fname_woext)

    #pd.set_option('display.max_columns', 500)
    r,n = make_coopfile(infile)
    r.to_csv("coop_%s.tsv" % fname_woext,index=False,sep="\t")
    n.to_csv("negctrl_%s.tsv" % fname_woext,index=False,sep="\t")

    print("Making training file...")
    """
    # Making the classification file
    r = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/coop_array_files/coop_20191004_258614510001_ETS1_550_5_1-4_alldata.tsv",sep="\t") ###
    #pd.set_option('display.max_columns', None)
    n = pd.read_csv("/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/coop_array_files/negctrl_20191004_258614510001_ETS1_550_5_1-4_alldata.tsv",sep="\t") ###
    """
    probedata = ProbeData(r,n,percentile=negcutoff)

    classification = classifier.classify_per_orientation(probedata, pvalthres)
    utils.dictlist2file(classification,"cooplabeled_%s.txt" % fname_woext,listval=True)
    #classification = utils.read_dictlist_file("/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191030_coop-PBM_Ets1_v1_2nd/coop_array_files/cooplabeled_20191004_258614510001_ETS1_550_5_1-4_alldata.txt", as_int=True)
    #idx_list = set(classification['cooperative_o1'] + classification['cooperative_o2'] +  classification['additive_o1'] +  classification['additive_o2']+  classification['anticoop_o1'] +  classification['anticoop_o2'])
    #idx_list = random.sample(idx_list,1000)
    #idx_list.sort()
    #probedata.scatter_boxplot_row(5976,log=True)
    #probedata.multi_scatter_boxplot(idx_list,log=True)

    print("Make scatter plot for each inconsistent classification")
    class_main = ["additive_o1","additive_o2","cooperative_o1","cooperative_o2","anticoop_o1","anticoop_o2"]
    subset = {k:classification[k] for k in classification if k not in class_main}
    """
    for sub in subset:
        print("    %s" % sub)
        probedata.multi_scatter_boxplot(subset[sub],log=True,filepath="boxplot-%s.pdf" % sub)
    """

    print("Plot median binding sum for all orientations...")

    cmain = {k:classification[k] for k in classification if k in class_main + ["cooperative_overlap","additive_overlap","anticoop_overlap"]}
    classifier.plot_median_binding_sum(probedata,cmain,1,log=True,tfname=tfname,plotname="all_plot_o1_log_%d" % negcutoff) # ,plotnonsignif=False
    classifier.plot_median_binding_sum(probedata,cmain,2,log=True,tfname=tfname,plotname="all_plot_o2_log_%d" % negcutoff)
    classifier.plot_median_binding_sum(probedata,cmain,0,log=True,tfname=tfname,plotname="all_plot_o1_2_log_%d" % negcutoff)


    print("Make training data") # need to make df_all
    # We only need classification regardless of the inconsistency here
    classification_main = {k:classification[k] for k in classification if k in class_main}
    class_main1 = ["cooperative_overlap","additive_overlap","anticoop_overlap"]
    classification_main1 = {k:classification[k] for k in classification if k in class_main1}
    class_main2 = ["cooperative_o1_anticoop_o2","cooperative_o2_anticoop_o1"]
    classification_main2 = {k:classification[k] for k in classification if k in class_main2}

    df_all = make_training(probe_analysis_path, "mutated_probes.*\.(tsv|csv)$", probedata, classification_main)
    df_all.to_csv("training_all_%s.tsv" % fname_woext, sep="\t", index=False)
    df_overlap1 = make_training(probe_analysis_path, "mutated_probes.*\.(tsv|csv)$", probedata, classification_main1)
    df_overlap1.to_csv("training_overlap_%s.tsv" % fname_woext, sep="\t", index=False)
    df_overlap2 = make_training(probe_analysis_path, "mutated_probes.*\.(tsv|csv)$", probedata, classification_main2)
    df_overlap2.to_csv("training_with_coop_anti_%s.tsv" % fname_woext, sep="\t", index=False)
    #pd.set_option('display.max_columns', 500)
    #df.loc[(df['label'] == 'cooperative') | (df['label'] == 'anticoop')].to_csv("training_cooperative_anticoop.tsv", sep="\t", index=False)
    #df.loc[(df['label'] == 'cooperative') | (df['label'] == 'additive')].to_csv("training_cooperative_additive.tsv", sep="\t", index=False)

    """
    # Making the plots for original genomic sequences
    print("----- Labeling and plotting of the original genomic sequences -----")
    # 24824
    r_orig = r[~r['Name'].str.contains(r'dist|weak')] # ets1_A549_seq26_dist_5', 'ets1_A549_seq26_dist_5_dup1'
    probedata_orig = ProbeData(r_orig,n,percentile=negcutoff)
    classification_orig = classifier.classify_per_orientation(probedata_orig, pvalthres, classify_inconsistency=False)
    classifier.plot_median_binding_sum(probedata_orig,classification_orig,1,log=True,tfname=tfname,plotname="orig_plot_o1_log_%d" % negcutoff)
    classifier.plot_median_binding_sum(probedata_orig,classification_orig,2,log=True,tfname=tfname,plotname="orig_plot_o2_log_%d" % negcutoff)
    classifier.plot_median_binding_sum(probedata_orig,classification_orig,0,log=True,tfname=tfname,plotname="orig_plot_o1_2_log_%d" % negcutoff)

    # Making the plots for custom sequences
    print("----- Labeling and plotting of the custom sequences -----")
    r_cust = r[r['Name'].str.contains(r'dist|weak')]
    probedata_cust = ProbeData(r_cust,n,percentile=negcutoff)
    classification_cust = classifier.classify_per_orientation(probedata_cust, pvalthres, classify_inconsistency=False)
    classifier.plot_median_binding_sum(probedata_cust,classification_cust,1,log=True,tfname=tfname,plotname="cust_plot_o1_log_%d" % negcutoff) # ,plotnonsignif=False
    classifier.plot_median_binding_sum(probedata_cust,classification_cust,2,log=True,tfname=tfname,plotname="cust_plot_o2_log_%d" % negcutoff)
    classifier.plot_median_binding_sum(probedata_cust,classification_cust,0,log=True,tfname=tfname,plotname="cust_plot_o1_2_log_%d" % negcutoff)
    print("Done!")
    """
