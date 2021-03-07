import pandas as pd
import chip2probe.modeler.dnashape as ds
from chip2probe.util import bio
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def align(seq, move):
    """
    move, if negative, move to the left (append right)
    """
    mv = abs(move)
    if move > 0:
        return bio.gen_randseq(mv) + seq[:-mv]
    else:
        return seq[mv:] + bio.gen_randseq(mv)

if __name__ == "__main__":
    params = {'axes.labelsize': 22,
          'axes.titlesize': 18,
          "xtick.labelsize" : 15, "ytick.labelsize" : 15 , "axes.labelsize" : 14}
    plt.rcParams.update(params)
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    df = pd.read_csv("%s/ch1_ch2/training_pwm_ori.tsv" % basepath, sep="\t")
    dist = sorted(df["distance"].unique().tolist())[1:]
    oris = df["orientation"].unique().tolist()

    fastadict = dict(zip(df["Name"], df["sequence"]))
    ds = ds.DNAShape(fastadict)

    #df[(df['distance'] == 5) & (df['orientation'] == "-/-")].to_csv("seq.csv",index=False)
    #df = df[df['distance'] == 5].sort_values("orientation").to_csv("seqtest.csv",index=False)
    print(df.groupby(["orientation","label"])["Name"].count())
    n = 0
    for ori in ["-/-"]:
        #curdf = pd.DataFrame(df[df['distance'] == d])
        curdf = pd.DataFrame(df[(df['orientation'] == ori)])
        s1 = int(curdf['ets_pos'].mode())
        s2 = curdf[curdf['ets_pos'] == s1].iloc[0]['runx_pos']
        curdf["seqalign"] = curdf.apply(lambda x: align(x['sequence'], s1-x['ets_pos']), axis=1)
        curdf = curdf[curdf["label"] == "cooperative"]
        ss = ""
        for s in curdf["seqalign"].tolist():
            ss += s + "\n"
        with open("ori.txt",'w') as f:
            f.write(ss)
    with PdfPages("save.pdf") as pdf:
        for ori in oris:
            #curdf = pd.DataFrame(df[df['distance'] == d])
            curdf = pd.DataFrame(df[df['orientation'] == ori])
            s1 = int(curdf['ets_pos'].mode())
            s2 = curdf[curdf['ets_pos'] == s1].iloc[0]['runx_pos']
            curdf["seqalign"] = curdf.apply(lambda x: align(x['sequence'], s1-x['ets_pos']), axis=1)
            fig = plt.figure(figsize=(12,12))
            print(curdf)
            ds.plot_average(curdf, linemarks=[s1, s2], base_count=False, in_fig=fig, lblnames=["cooperative","independent"], pltlabel="orientation %s"%ori, label_title=True)
            pdf.savefig(fig)
            plt.close()
