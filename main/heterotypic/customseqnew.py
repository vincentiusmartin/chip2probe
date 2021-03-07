from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import chip2probe.training_gen.arranalysis as arr

import chip2probe.util.stats_r as st

def correct_type(type, ori, orig_ori):
    flip = {"reduce_dist_from_glass":"reduce_dist_from_free", "reduce_dist_from_free":"reduce_dist_from_glass"}
    if type != "original" and ori != orig_ori:
        return flip[type]
    else:
        return type

def read_filter(path, key, info, seqdict, coltype, sep="\t"):
    df = pd.read_csv(path, sep=sep) \
           .merge(info, on="Name")
    df = df[df["Name"].str.contains(key,na=False) | (df["type"] == "original")]
    df["seqid"] = df["ori_seq"].map(seqdict)
    df["Name"] = df.apply(lambda x: "seq%d-d%d" % (x["seqid"], x[coltype]), axis = 1)
    asc = True if coltype == "step" else False
    df["type"] = df.apply(lambda x: correct_type(x["type"], x["ori"], x["orig_ori"]), axis=1)
    df = df.sort_values(["seqid", coltype], ascending=[True,asc])
    return df.drop_duplicates()

def print_shifted_pred(df, coltype):
    asc = True if coltype == "step" else False
    df_sorted = df.sort_values(by=["seqid",coltype], ascending=[True,asc])
    types = set(lbls["type"].unique()) - {'original'}
    for t in types:
        change_dict = {"er":{}, "re":{}}
        df_cur = df_sorted[df_sorted["type"] == t]
        grouped = df_cur.groupby('seqid')
        for name, group in grouped:
            first = group.head(1).iloc[0]
            last = group.tail(1).iloc[0]
            for o in ["er", "re"]:
                if (first[o] == "cooperative" or first[o] == "additive") and \
                    (last[o] == "cooperative" or last[o] == "additive") and \
                    (first[o] != last[o]):
                    change_str = "%s,%s" % (first[o], last[o])
                    if change_str not in change_dict[o]:
                        change_dict[o][change_str] = 1
                    else:
                        change_dict[o][change_str] += 1
        print(t, change_dict)

def plot_label_change(lbls, path, types, coltype):
    asc = True if coltype == "step" else False
    label_map = {'below_cutoff':-2, 'anticoop':-1, 'additive':0,  'ambiguous':1, 'cooperative':2}
    idx_int, idx_str =  [-2,-1,0,1,2], ['<_cut', 'anticoop', 'additive',  'ambig', 'coop']
    n = 0
    lbls = lbls.sort_values("seqid")
    seqids = list(lbls["seqid"].unique())

    with PdfPages(path) as pdf:
        for s in seqids:
            for t in types:
                l = lbls[ (lbls["seqid"] == s) &
                       ((lbls["type"] == "original") | (lbls["type"] == t))
                    ].reset_index(drop=True).sort_values(by=coltype, ascending=asc)
                if l.empty or l.shape[0] == 1:
                    continue
                if n == 0:
                    fig = plt.figure(figsize=(18,18))
                    fig.subplots_adjust(hspace=0.4,wspace=0.5)
                n += 1
                ax = fig.add_subplot(4,3,n)
                l["er"] = l['er'].map(label_map).astype('int64')
                l["re"] = l['re'].map(label_map).astype('int64')
                l.plot(x=coltype, y=["er","re"], ax=ax, marker='o', color=["red", "blue"])
                d = list(l[coltype])
                ax.set_xticks(d)
                ax.set_xlim(d[0],d[-1])
                ax.set_yticks(idx_int)
                ax.set_yticklabels(idx_str)
                ax.set_title("sequence %d - %s" % (s,t), fontsize=8)
                if n == 4*3:
                    pdf.savefig(fig)
                    plt.close()
                    n = 0
        pdf.savefig(fig)
        plt.close()
    return 0

def plot_confidence_change(indiv, two, coltype):
    asc = True if coltype == "step" else False
    indivsub = indiv[["seqid","type",coltype,"ori","affinity"]]
    twosub = two[["seqid","type",coltype,"ori","affinity"]]
    indiv_g = indivsub.groupby(["seqid","type",coltype, "ori"])
    two_g = twosub.groupby(["seqid","type",coltype, "ori"])
    df_p = two_g.apply(lambda g: pd.Series({"p":st.wilcox(g["affinity"].tolist(), indiv_g.get_group(g.name)["affinity"].tolist(), "greater")})).reset_index()
    df_p["p"] = np.log(df_p["p"])
    df_med = two_g.apply(lambda g: pd.Series({"med_delta":g["affinity"].median() - indiv_g.get_group(g.name)["affinity"].median()})).reset_index()
    df = df_p.merge(df_med, on=["seqid",coltype,"ori","type"]).sort_values(["seqid","distance"], ascending=[True,asc])
    types = list(set(indiv["type"].unique()) - {'original'})
    dfg = df.groupby("seqid")
    n = 0
    for toplot in ["med_delta", "p"]:
        with PdfPages("change_%s.pdf" % toplot) as pdf:
            dfg = df.groupby("seqid")
            for name, g in dfg:
                for t in types:
                    if n == 0:
                        fig = plt.figure(figsize=(18,18))
                        fig.subplots_adjust(hspace=0.4,wspace=0.5)
                    n += 1
                    ax = fig.add_subplot(4,3,n)
                    cur_g = g[(g["type"] == t) | (g["type"] == "original")]
                    if cur_g.empty or cur_g.shape[0] == 1:
                        continue
                    d = list(cur_g["distance"].unique())
                    cur_g.set_index('distance', inplace=True)
                    cur_g.groupby('ori')[toplot].plot(legend=True, ax=ax)
                    ax.set_xticks(d)
                    ax.set_xlim(d[0],d[-1])
                    ax.set_title("sequence %d - %s" % (name,t), fontsize=8)
                    if n == 4*3:
                        pdf.savefig(fig)
                        plt.close()
                        n = 0
            pdf.savefig(fig)
            plt.clf()
            plt.close()

def plot_each_probe(indiv, two, coltype): #lbls
    print("Plotting each probe")
    types = set(indiv["type"].unique()) - {'original'}
    #plot_label_change(lbls, "lblchange.pdf", types, coltype)

    for t in types:
        cur_indiv = pd.DataFrame(indiv[(indiv["type"] == t) | (indiv["type"] == 'original')])
        cur_two = two[(two["type"] == t) | (two["type"] == 'original')]
        arr.plot_ori_inconsistency(cur_indiv, cur_two, log=True, prefix_path="%s_"%t, fixed_ax=True, namecol="Name",thresline=cutoff)

def plot_scatter_per_category(indiv, two, log=True):
    types = list(set(indiv["type"].unique()) - {'original'})
    dfmed = {}
    dfmed["indiv"] = indiv.groupby(["Name","ori","type"]).median().reset_index()[["Name","ori","type","affinity"]]
    dfmed["two"] = two.groupby(["Name","ori","type"]).median().reset_index()[["Name","ori","type","affinity"]]
    if log:
        dfmed["indiv"]["affinity"] = np.log(dfmed["indiv"]["affinity"])
        dfmed["two"]["affinity"] = np.log(dfmed["two"]["affinity"])
    oris = ["er","re"]
    for ch in ["indiv","two"]:
        print("R2 for %s" % ch)
        for o in oris:
            for t in types:
                df1 = dfmed[ch][(dfmed[ch]["type"] == t) & (dfmed[ch]["ori"] == o)][["Name","affinity"]]
                df2 = dfmed[ch][(dfmed[ch]["type"] == "shorter") & (dfmed[ch]["ori"] == o)][["Name","affinity"]]
                r2 = arr.plot_chamber_corr(df1, df2,  xlab="Intensity %s" % t, ylab="Intensity shorter",
                    valcol="affinity", title="%s, %s vs shorter, ori %s" % (ch,t,o), path="scatter_%s_%s_%s.png" % (ch,t,o))
                print("Orientation %s: %s vs shorter,  %.2f" % (o, t, r2))

def plot_delta_ch1ch2(indiv,two):
    indivmed = indiv[["seqid","ori","type","distance","affinity"]] \
        .groupby(["seqid","ori","type","distance"]) \
        .median().reset_index().sort_values(["seqid","distance"], ascending=[True,False])
    twomed = two[["seqid","ori","type","distance","affinity"]] \
        .groupby(["seqid","ori","type","distance"]) \
        .median().reset_index().sort_values(["seqid","distance"], ascending=[True,False])
    types = list(set(indivmed["type"].unique()) - {'original'})
    oris = list(indivmed["ori"].unique())
    df = indivmed.merge(twomed, on=["seqid","ori","type","distance"], suffixes=('_indiv', '_two'))
    df['diff'] = df["affinity_two"] - df["affinity_indiv"]
    for t in types:
        curdf = df[(df["type"] == t) | (df["type"] == "original")]
        for o in oris:
            curdf_ori = curdf[curdf["ori"] == o]
            #print("Number of sequences %s,%s: %d" % (t,o,len(curdf["seqid"].unique().tolist())))
            d = list(curdf_ori["distance"])
            y = list(curdf_ori["diff"])
            maxd_aff, mind_aff = curdf_ori[curdf_ori["distance"] == max(d)]['diff'].tolist(), curdf_ori[curdf_ori["distance"] == min(d)]['diff'].tolist()
            ax = plt.figure().add_subplot(111)
            plt.scatter(d, y)
            plt.xlim(max(d), min(d))
            slope, intercept, r, p, std_err = stats.linregress(d, y) #p-val with Wald Test
            plt.text(0.5, 0.5, 'p = %0.5f' % (p), color='red', fontsize=12, transform = ax.transAxes)
            print("%s %s ,p=%0.5f" % (t,o,p))
            dist_unique = list(set(d))
            abline_values = [slope * i + intercept for i in dist_unique]
            plt.title("%s_%s, slope %.2f" % (t,o,slope))
            plt.ylabel("ch2_intensity - ch1_intensity")
            plt.xlabel("distance from the second site start to the glass")
            plt.xticks(df["distance"].unique().tolist())
            plt.plot(dist_unique, abline_values, color='black')
            plt.tight_layout()
            plt.savefig("distance_site2_%s_%s.png" % (t,o))
            plt.clf()

def plot_delta_dist(indiv):
    indivsub = indiv[["seqid","ori","type","distance","affinity"]]
    indivmed = indivsub.groupby(["seqid","ori","type","distance"]).median().reset_index().sort_values(["seqid","distance"], ascending=[True,False])
    types = list(set(indiv["type"].unique()) - {'original'})
    oris = list(indivmed["ori"].unique())
    for t in types:
        curdf = indivmed[(indivmed["type"] == t) | (indivmed["type"] == "original")]
        for o in oris:
            curdf_ori = curdf[curdf["ori"] == o]
            print(curdf_ori)
            original_seqs = curdf[curdf["type"] == "original"][["seqid","ori","affinity"]].rename({"affinity":"diff"}, axis=1)
            curdf_ori = curdf_ori.merge(original_seqs, on=["seqid","ori"])
            curdf_ori["diff"] = curdf_ori["diff"] - curdf_ori["affinity"]
            d = list(curdf_ori["distance"])
            y = list(curdf_ori["diff"])
            plt.scatter(d, y)
            plt.xlim(max(d), min(d))
            slope, intercept = np.polyfit(d, y, 1)
            # Create a list of values in the best fit line
            dist_unique = list(set(d))
            abline_values = [slope * i + intercept for i in dist_unique]
            plt.title("%s_%s, slope %.2f" % (t,o,slope))
            plt.xticks(df["distance"].unique().tolist())
            plt.ylabel("intensity")
            plt.xlabel("distance")
            plt.plot(dist_unique, abline_values, color='black')
            plt.savefig("distance_effect_%s_%s.png" % (t,o))
            plt.clf()

def read_info(info_path):
    info = pd.read_csv(info_path)
    info["Name"] = info["Name"].apply(lambda x: "_".join(x.split("_")[:-3]))
    # info["distance"] = info['ets_start'] - info['step']
    info_orig = pd.DataFrame(info[info["type"] == "original"])
    #info_orig["distance"] = info["ets_start"] - info["step"]
    info_orig["distance"] = info.apply(lambda x: max(x["ets_start"], x["runx_start"]) - x["step"] ,axis=1)
    info_orig = info_orig[["ori_seq", "distance"]].drop_duplicates()
    info = info.merge(info_orig, on="ori_seq")
    info = info[["Name","distance","type","ori_seq","orig_ori"]]
    return info

if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/customseq"
    ch = "ch1_ch2"
    key = "_new_customseq"
    cutoff = np.log(860.98) # 446.035 np.log(860.98)
    pd.set_option("display.max_columns",None)
    coltype = "distance" # or step

    info = read_info("%s/custom_info_fz_worigori.csv" % basepath)
    oriseq_dict = {k: v for v, k in enumerate(info["ori_seq"].unique())}
    #info = info[info['ori_seq'] == "TTCTCTGCGGTGGCCAGAGTCAGAAGCGGATAAACA"] ###

    indiv = read_filter("%s/%s/indiv_df.tsv" % (basepath,ch),key, info, oriseq_dict, coltype)
    two = read_filter("%s/%s/two_df.tsv" % (basepath,ch),key, info, oriseq_dict, coltype)
    # indiv.to_csv("aaa.csv", index=False)
    # oriseqs = indiv["ori_seq"].unique()
    #
    plot_delta_ch1ch2(indiv,two)
    # plot_delta_dist(indiv)
    # print_shifted_pred(lbls, coltype)
    # plot_scatter_per_category(indiv, two)
    # plot_confidence_change(indiv, two, coltype)


    # 1. Plot each probe to see the changes
    # plot_each_probe(indiv, two, coltype)
