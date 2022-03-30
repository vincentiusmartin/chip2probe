import pandas as pd
import chip2probe.modeler.plotlib as plot
import seaborn as sns
import matplotlib.pyplot as plt
import chip2probe.util.stats as st
from decimal import Decimal
import numpy as np
import seaborn as sns

pd.set_option("display.max_columns", None)

maintf, cooptf = "runx1","ets1"
# er "#b22222","#FFA07A", re "#0343df","#75bbfd"
colors = ["#0343df","#75bbfd"]
df = pd.read_csv("output/Runx1Ets1/training/train_runx1_ets1.tsv", sep="\t")[["distance","orientation","%s_score" % maintf, "%s_score" % cooptf, "label"]]

# pd.set_option("display.max_rows",None)
# sz = df.groupby(["distance","orientation","label"]).size()
# print(sz)

maintf_lbl = '%s strength\n(main TF)' % maintf.capitalize()
cooptf_lbl = '%s strength\n(cooperator TF)' % cooptf.capitalize()
df.rename(columns={'%s_score' % maintf: maintf_lbl, '%s_score' % cooptf: cooptf_lbl}, inplace=True)

dfcoop = df[df["label"] == "cooperative"]
dfcoop = dfcoop.groupby(by=["distance","orientation"]).filter(lambda x: len(x) >= 5)
selected_dist_ori = dfcoop[["distance","orientation"]].drop_duplicates().sort_values(["distance","orientation"])
do_list = ["%d\n%s"%(x[0],x[1]) for x in selected_dist_ori.values.tolist()]
df = df.merge(selected_dist_ori, on=("distance","orientation")).sort_values(["distance","orientation"])
df = df[df["label"] == "cooperative"] ###

plt.rcParams["figure.figsize"] = (20,5.5)
df["xtick"] = df.apply(lambda x: "%d\n%s"%(x["distance"],x["orientation"]),axis=1)

col_y = '%s strength\n(cooperator TF)' % cooptf.capitalize()
grouped_coop = df[df["label"] == "cooperative"].groupby(by=["xtick"])[col_y].apply(list).to_dict()
grouped_indep = df[df["label"] == "independent"].groupby(by=["xtick"])[col_y].apply(list).to_dict()

"""
ax = sns.boxplot(data=df, x="xtick", y=col_y,
                hue="label", width=0.9, palette=colors)
ax.set_ylim(top=18)

pdict = {k:st.wilcox(grouped_coop[k],grouped_indep[k],alternative="greater") for k in do_list if k in grouped_coop and k in grouped_indep}
for i in range(len(do_list)):
    k = do_list[i]
    if not (k in grouped_coop and k in grouped_indep):
        continue
    p = st.wilcox(grouped_coop[k],grouped_indep[k],alternative="greater")
    y = 1.075*max(grouped_coop[k]), 1.075*max(grouped_indep[k])
    maxy = max(y[0],y[1])
    if p < 0.01:
        if p < 1e-10:
            annotate_text = "***"
            minus_pos = 0.4
        elif p < 1e-5:
            annotate_text = "**"
            minus_pos = 0.25
        else:
            annotate_text = "*"
            minus_pos = 0.12
        ax.annotate(annotate_text, xy=(i-minus_pos,maxy), zorder=0)
    pt1,pt2 = 1.05*max(grouped_coop[k]), 1.05*max(grouped_indep[k])
    plt.plot((i-0.25,i-0.25,i+0.22,i+0.22), (pt1,maxy,maxy,pt2), lw=1, c="#2680ff")


# df = df[df["label"] == "cooperative"]
# ax = sns.boxplot(data=df, x="xtick", y='%s strength\n(cooperator TF)' % cooptf.capitalize(),width=.6,color=colors[0])
# ax.set_ylim(-8,23)
# ax.set_yticks([-5, 0, 5, 10, 15])
# ax.get_legend().remove()
# plt.title("Cooperative Runx1 (cooperator TF) strength")
# plt.savefig("box.png")
"""

filtered_grp = {k: list(v) for k, v in df.groupby(["distance","orientation"])['%s strength\n(cooperator TF)' % cooptf.capitalize()]}

allps = {y: st.wilcox(filtered_grp[(5, '+/+')], filtered_grp[y], alternative="less") for y in list(filtered_grp.keys())}
# allps = {k: v for k, v in sorted(allps.items(), key=lambda item: item[1])}
for x in allps:
    pval = '{:0.2e}'.format(allps[x]) if allps[x] < 0.1 else "%.2f"%allps[x]
    print(x,":", pval)

"""
p1 = st.wilcox(cur_group["+/+"],cur_group["+/-"],alternative="less")
p2 = st.wilcox(cur_group["+/+"],cur_group["-/+"],alternative="less")
p3 = st.wilcox(cur_group["+/+"],cur_group["-/-"],alternative="less")
print("%.3e" % Decimal(p1),"%.3e" % Decimal(p2),"%.3e" % Decimal(p3))

g1 = df[(df["distance"] == 5) & (df["orientation"] == "+/+")]
g2 = df[~((df["distance"] == 5) & (df["orientation"] == "+/+"))]

for lbl in ["cooperative","independent"]:
    g1_c = g1[g1["label"] == lbl][maintf_lbl]
    g2_c = g2[g2["label"] == lbl][maintf_lbl]
    dat = [g1_c.tolist(),g2_c.tolist()]
    print(lbl, len(dat[0]), len(dat[1]))
    p_gr = st.wilcox(g1_c,g2_c,alternative="greater")
    p_ls = st.wilcox(g1_c,g2_c,alternative="less")
    p = p_gr if p_gr < p_ls else p_ls
    print("p",p)

    plt.rcParams["figure.figsize"] = (4,4)
    plt.yticks([-5, 0, 5, 10, 15])
    plt.ylim(-9,18)
    sns.boxplot(data=dat, width=.6, color=colors[1] if lbl == "independent" else colors[0])
    plt.savefig("fig_%s.pdf"%lbl)
    plt.clf()

# df = df[df["distance"] == 20]
# grp_ori = df.groupby("orientation")
# print(grp_ori["label"].value_counts())
# for idx,group in grp_ori:
#     name_id = idx.replace("/","")
#     plot.plot_box_categories(group, path="boxplot_%s.png" % name_id, incols=["%s strength\n(main TF)" % maintf.capitalize(), "%s strength\n(cooperator TF)" % cooptf.capitalize()], alternative="smaller", color=color)
"""
