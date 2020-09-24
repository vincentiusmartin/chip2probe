
import pandas as pd

if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/ch1_ch2"
    chdf = pd.read_csv("%s/coop_ch1_vs_ch2/sequence_labeled_p06.tsv" % basepath, sep="\t")[["Name","label"]].drop_duplicates()
    mtdf = pd.read_csv("%s/coop_wt_vs_mt/EtsRunx_wt_cutoff/name_labeled.csv" % basepath)
    mtdf = mtdf[mtdf["label"] != "below_cutoff"]

    print("Overlap in wtmt and ch1ch2")
    probe_overlap = chdf.merge(mtdf, on="Name")
    match_lbl = probe_overlap[probe_overlap["label_x"] == probe_overlap["label_y"]]
    print(probe_overlap.shape[0])
    print(match_lbl["label_x"].value_counts())

    print("Only in ch1ch2")
    ch_only_names = set(chdf["Name"]).difference(mtdf["Name"])
    ch_only_df = chdf[chdf["Name"].isin(ch_only_names)]
    print(ch_only_df["label"].value_counts())

    print("Only in wtmt")
    mt_only_names = set(mtdf["Name"]).difference(chdf["Name"])
    mt_only_df = mtdf[mtdf["Name"].isin(mt_only_names)]
    print(mt_only_df["label"].value_counts())
