import pandas as pd
import chip2probe.modeler.plotlib as plot
import chip2probe.training_gen.arranalysis as arr

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

"""
1. original sequence
2. fix site2, remove mid, bring site 1 closer to site2, adding flank on the left
3. fix site1, remove mid, bring site 2 closer to site 1, adding flank on the right
4. remove mid, bring both sites, no flank added
"""

def read_filter(path, key, sep="\t"):
    df = pd.read_csv(path, sep=sep)
    df = df[df["Name"].str.contains(key,na=False)]
    df["Name"] = df["Name"].str.replace(key,"")
    df["idx"] = df["Name"].str.split("_").str[0].str[3:].astype('int64')
    df = df.sort_values(by=["idx"])
    return df

def read_combine_filter(cust, ref, key):
    cust_df = pd.read_csv(cust)
    ref_df = pd.read_csv(ref)[["Name", "ori_seq","step","type"]]
    ref_df["Name"] = ref_df["Name"].apply(lambda x: "_".join(x.split("_")[:-3]))
    cmb = cust_df.merge(ref_df, on="Name")
    cmb = cmb[(cmb["step"] == 0) | cmb["Name"].str.contains(key,na=False)].sort_values(by=["ori_seq"])
    cmb["Name"] = cmb["Name"].str.replace("_clean_customseq","").str.replace("_new_customseq","")
    cmb["idx"] = cmb["Name"].str.split("_").str[0].str[3:].astype('int64')
    cmb = cmb.sort_values(by=["idx"])
    cmb.to_csv("lalala.csv", index=False)

def assign_type(x):
    if x["step"] == 0:
        return "orig"
    elif x["step"] > 0 and len(x["Sequence"]) == 36:
        return "to_right"
    elif x["step"] < 0 and len(x["Sequence"]) == 36:
        return "to_left"
    else:
        return "shorter"

def plot_label_change(lbls,titles,path):
    label_map = {'below_cutoff':-2, 'anticoop':-1, 'additive':0,  'ambiguous':1, 'cooperative':2}
    idx_int, idx_str =  [-2,-1,0,1,2], ['<_cut', 'anticoop', 'additive',  'ambig', 'coop']
    n = 0
    with PdfPages(path) as pdf:
        for i in range (0,len(lbls)):
            l = lbls[i].set_index("step")
            if n == 0:
                fig = plt.figure(figsize=(18,18))
                fig.subplots_adjust(hspace=0.4,wspace=0.5)
            n+=1
            ax = fig.add_subplot(4,4,n)
            l["er"] = l['er'].map(label_map).astype('int64')
            l["re"] = l['re'].map(label_map).astype('int64')
            l[["er","re"]].plot.line(ax=ax)
            #ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.set_xticks(l.index)
            ax.set_yticks(idx_int)
            ax.set_yticklabels(idx_str)
            ax.set_title(titles[i], fontsize=8)
            if n == 4*4:
                pdf.savefig(fig)
                plt.close()
                n = 0
        pdf.savefig(fig)
        plt.close()
    return 0

if __name__=="__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/customseq"
    ch = "ch3_ch4"

    # x = pd.read_csv("%s/%s/df_custom_ch1.csv" % (basepath,ch))
    # x = x[x["Name"].str.contains("custom",na=False)]
    # x["idx"] = x["Name"].str.split("_").str[0].str[3:].astype('int64')
    # x.sort_values(by="idx").to_csv("custom_raw.csv",index=False)

    key = "_clean_customseq"
    seqstep = pd.read_csv("%s/clean_custom.csv" % basepath)[["sequence", "step"]].rename({"sequence":"Sequence"}, axis=1).drop_duplicates()
    custom_df = read_filter("%s/%s/df_custom_ch3.csv" % (basepath,ch), key, ",")
    #read_combine_filter("%s/%s/df_custom_ch1.csv" % (basepath,ch), "%s/custom_all.csv" % basepath, key)

    name_step_map = custom_df[["Name","Sequence"]].merge(seqstep, on="Sequence", how="left").dropna()
    name_step_map["step"] = name_step_map["step"].astype('int64')
    name_step_map["type"] = name_step_map.apply(lambda x: assign_type(x), axis=1)
    name_step_map = name_step_map[["Name","step","type"]].drop_duplicates()
    custom_df = custom_df.merge(name_step_map, on="Name").sort_values(by=["idx"])

    indiv_df = read_filter("%s/%s/indiv_df.tsv" % (basepath,ch), key).merge(name_step_map, on="Name")
    two_df = read_filter("%s/%s/two_df.tsv" % (basepath,ch), key).merge(name_step_map, on="Name")

    lbldf = read_filter("%s/%s/lbl_df.tsv" % (basepath,ch), key).merge(name_step_map, on="Name")
    #lbldf["Name"] = lbldf.apply(lambda x: "%s_step%d" % (x["Name"], x["step"]), axis=1)
    orig_lbl = lbldf[lbldf["step"] == 0]
    l = list(orig_lbl.index) + [lbldf.shape[0]]
    seqidlist =  [y for x in [[n]*(l[n+1]-l[n]) for n in range(len(l)-1)] for y in x]
    unique_id = [n for n in range(len(l)-1)]
    lbldf["seqid"] = seqidlist
    name_id_map = lbldf[["Name","seqid","step", "type"]]
    name_id_map["fullname"] = name_id_map.apply(lambda x: "%s_step%d" % (x["seqid"], x["step"]), axis=1)
    name_id_map = name_id_map[["Name","fullname"]]
    indiv_df = indiv_df.merge(name_id_map, on="Name")
    two_df = two_df.merge(name_id_map, on="Name")

    for t in ["to_left", "to_right", "shorter"]:
        lbl = lbldf[(lbldf["type"] == "orig") | (lbldf["type"] == t)].reset_index(drop=True)
        orig_lbl_sub = lbl[lbl["step"] == 0]
        lsub = list(orig_lbl_sub.index) + [lbl.shape[0]]
        list_lbls = [lbl.iloc[lsub[n]:lsub[n+1]] for n in range(len(lsub)-1)]
        plot_label_change(list_lbls, unique_id, "lblchange_%s.pdf" % t)

        ind = indiv_df[(indiv_df["type"] == "orig") | (indiv_df["type"] == t)]
        two = two_df[(two_df["type"] == "orig") | (two_df["type"] == t)]
        # ind = ind.reindex(ind["step"].abs().sort_values().index)
        # two = two.reindex(two["step"].abs().sort_values().index)
        # ind["Name"] = ind["step"]
        # two["Name"] = two["step"]
        arr.plot_ori_inconsistency(ind, two, log=True, prefix_path="%s_"%t, fixed_ax=True, namecol="fullname")

    # both_ori_labeled = read_filter("%s/ch1_ch2/sequence_labeled_normalized.tsv" % basepath, key ,sep="\t")
    # custom_df = custom_df[custom_df["Name"].str.contains("new_customseq",na=False)].sort_values(by=["Name"])
    # custom_df["seqid"] = custom_df["Name"].str.split('_').str[0]
    # train_df = pd.read_csv("%s/ch1_ch2/training_pwm.tsv" % basepath, sep="\t")
    # train_df["seqid"] = train_df["Name"].str.split('_').str[0]
    #
    # custom_df = custom_df[["seqid","label"]].drop_duplicates()
    # orig_df = train_df[["seqid","label"]].drop_duplicates()
    # df_merge = custom_df.merge(orig_df, on="seqid",suffixes=("_cust", "_orig"))
    #
    # print(df_merge.groupby(["label_orig","label_cust"]).count())
    #
    # select = df_merge[(df_merge["label_orig"] == "cooperative") & (df_merge["label_cust"] == "independent")][["seqid"]].merge(train_df)
    # print(select.columns)
    #
    # plot.plot_box_categories(select, incols=["ets_score", "runx_score"], alternative="smaller")
