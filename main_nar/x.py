import pandas as pd

pd.set_option("display.max_columns",None)
# df = pd.read_csv("output/Ets1Runx1/label_pr/seqlbled_ets1_runx1.tsv", sep="\t")
# print(df.shape[0])
# print(df[df["label"] == "cooperative"])

# df = pd.read_csv("output/Ets1Ets1_v2/label_pr/wtm1m2m3_trimmed.csv")
# print(df.shape[0])
# print(df["label"].value_counts())

df = pd.read_csv("output/Ets1Runx1/training/train_ets1_runx1.tsv",sep="\t")
print(df.shape[0])
dfx = df[df["distance"] == 5]
print(dfx["label"].value_counts())
