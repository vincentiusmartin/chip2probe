import pandas as pd
pd.set_option("display.max_columns",None)

basepath = "output/Ets1Ets1_v2"
df = pd.read_csv("%s/training/train_ets1_ets1.tsv" % basepath, sep="\t")[["distance","orientation","site_str_score","site_wk_score","label"]]
df_count = pd.DataFrame(df.value_counts(subset=["distance","orientation","label"])).sort_values(by=["distance","orientation","label"]).reset_index()
df_count.columns = ["distance","orientation","label","count"]
df_count['Percentage'] = 100 * df_count['count'] / df_count.groupby(["distance","orientation"])['count'].transform('sum')

print(df_count)
for d in range(4,25):
    df_filter_dist = df_count[(df_count["distance"] == d)]
    print("d=%d" % d,"%.2f" % float(df_filter_dist[df_filter_dist["label"] == "cooperative"][["count"]].sum() / df_filter_dist[["count"]].sum()))

dft_med = df.groupby(["distance","orientation","label"])[['site_str_score','site_wk_score']].agg('median').reset_index()

df_conf = df_count.merge(dft_med, on=["distance","orientation","label"])
df_conf[df_conf["label"] == "cooperative"].to_csv("ets1_ets1_conf_coop.csv",float_format='%.2f')
df_conf[df_conf["label"] == "independent"].to_csv("ets1_ets1_conf_independent.csv",float_format='%.2f')
