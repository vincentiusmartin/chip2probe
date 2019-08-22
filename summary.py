import os, re
import pandas as pd


result_path = "../result"

pattern = "mutated_probes.*\.(tsv)$"

def useful_commands(df):
    # 1. Group by and count based on span + nsites
    df.groupby(["span","nsites","ecutoff"]).size().reset_index(name='counts')
    # 2. Get rows with duplicated wt sequence
    df[df.duplicated(['wt'])]
    # 3a. SELECT WHERE WT == ""
    df.loc[df['wt'] == ""]
    # 3b. Only get some columns
    df.loc[df['wt'] == ""][["nsites","span","tf","seqlabel"]]


def coop_probes_to_df(result_path, pattern):
    re_pattern = re.compile(pattern)
    agg = []
    for path, dirs, filenames in os.walk(result_path):
        for fname in filter(lambda name:re_pattern.match(name),filenames):
            tf_name = path.split("/")[-2] # should be the second from the last
            fsplit = os.path.splitext(fname)[0].split("_")
            span = fsplit[-1][len("span"):]
            nsites = fsplit[-2][1:]
            filepath = os.path.join(path, fname)

            probe_df = pd.read_csv(filepath, sep="\t")
            probe_df["tf"] = tf_name
            probe_df["span"] = span
            probe_df["nsites"] = nsites
            agg.extend(probe_df.to_dict('records'))
    df = pd.DataFrame(agg, columns=["wt","m1","m2","m3","ecutoff","nsites","egapthres","span","tf","key"], index=None)
    df = df.rename(columns={'key': 'seqlabel'})
    return df

def remove_wt_duplicates(df):
    newdf = df.drop_duplicates('wt', inplace = False)
    return newdf

if __name__=="__main__":
    summary = []
    df = coop_probes_to_df(result_path,pattern)
    no_dup = remove_wt_duplicates(df)
    no_dup.to_csv("out.csv")
    print(no_dup.groupby(["span","nsites","ecutoff","egapthres"]).size().reset_index(name='counts'))
    #print(no_dup.groupby(["nsites","span","ecutoff"]).size().reset_index(name='counts').to_string(index=False))
    #print(df.loc[df['tf'] == "ets1_HepG2"].groupby(["span","nsites","ecutoff","tf"]).size().reset_index(name='counts'))

    #df = pd.read_csv("../result/ets1_HepG2/analysis_result/mutated_probes_sites_within_d2_span50.tsv")
    print(df.shape,no_dup.shape)
