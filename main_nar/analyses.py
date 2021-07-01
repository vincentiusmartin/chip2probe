import pandas as pd
import chip2probe.modeler.plotlib as plot
import seaborn as sns
import matplotlib.pyplot as plt
import chip2probe.util.stats as st
from decimal import Decimal

pd.set_option("display.max_columns", None)

maintf, cooptf = "ets1","runx1"
color = ["#b22222","#FFA07A"]
df = pd.read_csv("output/Ets1Runx1/training/train_ets1_runx1.tsv", sep="\t")
# pd.set_option("display.max_rows",None)
sz = df.groupby(["distance","orientation","label"]).size()
# print(sz)

df.rename(columns={'%s_score' % maintf: '%s strength\n(main TF)' % maintf.capitalize(), '%s_score' % cooptf: '%s strength\n(cooperator TF)' % cooptf.capitalize()}, inplace=True)

"""

df = df[df["label"] == "cooperative"]


plt.rcParams["figure.figsize"] = (10,3)
grouped = df.groupby(by="distance")
cur_group = {elm[0]:list(elm[1]) for elm in grouped["Runx1 strength\n(cooperator TF)"]}
labels, data = [*zip(*cur_group.items())]
ax = sns.boxplot(data=data, width=.6, color="#b22222")
ax.set_xticklabels(labels)
plt.title("Cooperative Runx1 (cooperator TF) strength")
plt.savefig("box.png")

allps = {y: st.wilcox(cur_group[5], cur_group[y], alternative="less") for y in list(cur_group.keys())}
for x in allps:
    print(x,":", allps[x])


p1 = st.wilcox(cur_group["+/+"],cur_group["+/-"],alternative="less")
p2 = st.wilcox(cur_group["+/+"],cur_group["-/+"],alternative="less")
p3 = st.wilcox(cur_group["+/+"],cur_group["-/-"],alternative="less")
print("%.3e" % Decimal(p1),"%.3e" % Decimal(p2),"%.3e" % Decimal(p3))
"""

df = df[(df["distance"] == 5)]
grp_ori = df.groupby("orientation")
print(grp_ori["label"].value_counts())
for idx,group in grp_ori:
    name_id = idx.replace("/","")
    plot.plot_box_categories(group, path="boxplot_%s.png" % name_id, incols=["%s strength\n(main TF)" % maintf.capitalize(), "%s strength\n(cooperator TF)" % cooptf.capitalize()], alternative="smaller", color=color)
