import argparse
import pandas as pd
import chip2probe.training_gen.arranalysis as arr
from pathlib import Path
import chip2probe.util.bio as bio
from difflib import SequenceMatcher

def fix_naming(df):
    df = df.drop(columns=["ori"])
    df_ori = df[["Name","Sequence","type"]] \
        .drop_duplicates()
    df_ori["ori"] = df_ori \
        .groupby(["Name","type"],as_index=False) \
        .cumcount() + 1
    df_ori["ori"] = "o" + df_ori["ori"].astype(str)
    df = df.merge(df_ori, on=["Name","Sequence","type"])
    return df

def fixdup_single_name(nm, df):
    dseq1 = df[(df["type"] == "wt")][["Sequence"]].rename({"Sequence":"Sequence_x"}).drop_duplicates()
    dseq2 = df[["Sequence"]].rename({"Sequence":"Sequence_y"}).drop_duplicates()

    dseq1_o1_seqs = df[(df["type"] == "wt") & (df["ori"] == "o1")]["Sequence"].drop_duplicates().tolist()
    rcgrp = {}
    for i in range(len(dseq1_o1_seqs)):
        rcgrp[dseq1_o1_seqs[i]] = i
        rcgrp[bio.revcompstr(dseq1_o1_seqs[i])] = i

    dseq1["key"], dseq2["key"] = 1,1
    comb = dseq1.merge(dseq2, on="key").drop("key", 1)  # get cartesian product between the two

    comb["simscore"] = comb.apply(lambda row: SequenceMatcher(None, row["Sequence_x"], row["Sequence_y"]).ratio(), axis=1)
    comb = comb.sort_values(["Sequence_x","simscore"], ascending=[True,False])
    selectedcmb = comb.groupby("Sequence_x").head(4).rename(columns={"Sequence_y":"Sequence"})
    selectedcmb["group"] = selectedcmb["Sequence_x"].apply(lambda x: "%s_%s" % (nm, rcgrp[x]))
    selectedcmb = selectedcmb[["Sequence", "group"]].drop_duplicates()

    named = df.merge(selectedcmb, on="Sequence")
    return named

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Generate clean probe file with just the necessary information')
    parser.add_argument(action="store", dest="path", type=str, help='File input path')
    parser.add_argument('-k', '--key', action="store", dest="key", type=str,  help='include key string')
    parser.add_argument('-e', '--exkey', action="store", dest="exkey", type=str,  help='exclude key string')
    parser.add_argument('-n', '--negkey', action="store", dest="negkey", default="NegativeCtrl", type=str, help='include negctrl key string')
    parser.add_argument('-c', '--seqcols', action="store", dest="sc", default="Name,type,rep,ori", type=str, help='sequence cols order, str separated by comma')
    parser.add_argument('-d', '--negcols', action="store", dest="nc", default="Name,rep,ori", type=str, help='negctrl cols order, str separated by comma')
    parser.add_argument('-f', '--fixnaming', action="store_true", dest="fixnaming", help='fix naming error in the orientatiion files')
    parser.add_argument('-g', '--fixdup', action="store_true", dest="fixdup", help='fix naming duplicates')
    args = parser.parse_args()

    # python3 clean_file.py input/original/Ets1_Ets1.txt -k "ets1" -e "dist|weak" -g
    # python3 clean_file.py input/original/Ets1_only.txt -k "all_clean_seqs" -n "negative_controls" -f
    # python3 clean_file.py input/original/Ets1_Runx1.txt -k "all_clean_seqs" -n "negative_controls" -f
    # python3 clean_file.py input/original/Runx1_only.txt -k "all_clean_seqs" -n "negative_controls" -f
    # python3 clean_file.py input/original/Runx1_Ets1.txt -k "all_clean_seqs" -n "negative_controls" -f

    filename = Path(args.path).stem
    df, neg = arr.read_chamber_file(args.path, includekey=args.key, negkey=args.negkey, excludekey=args.exkey,
                seqcols=args.sc.split(","), negcols=args.nc.split(","))
    df.rename(columns={"Alexa488Adjusted": "intensity"}, inplace=True)
    neg.rename(columns={"Alexa488Adjusted": "intensity"}, inplace=True)

    if args.fixnaming:
        df = fix_naming(df)
    if args.fixdup:
        pd.set_option("display.max_columns",None)
        g = df.groupby("Name")
        dflist = []
        for nm, gdf in g:
            if gdf.shape[0] != 24:
                gdf = fixdup_single_name(nm, gdf)
                gdf["Name"] = gdf["group"]
                gdf = gdf.drop(columns=["group"])
            dflist.extend(gdf.to_dict('records'))
        df = pd.DataFrame(dflist)

    df.to_csv("%s_pr_clean.csv" % filename, index=False)
    neg.to_csv("%s_neg_clean.csv" % filename, index=False)
