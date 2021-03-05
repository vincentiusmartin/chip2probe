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
    df = pd.read_csv("output/Ets1Runx1/training/train_ets1_runx1.tsv", sep="\t")
    dist = sorted(df["distance"].unique().tolist())[1:]
    oris = df["orientation"].unique().tolist()

    fastadict = dict(zip(df["Name"], df["Sequence"]))
    ds = ds.DNAShape(fastadict)
    n = 0
    labels = ["independent", "cooperative"]

    with PdfPages("save.pdf") as pdf:
        for d in [5]:
            curdf = pd.DataFrame(df[df['distance'] == d])
            print(curdf.groupby("orientation").count()[["Name"]])
            for ori in ["+/+"]: # oris
                curdf = curdf[curdf['orientation'] == ori]
                s1 = int(curdf['ets1_pos'].mode())
                s2 = curdf[curdf['ets1_pos'] == s1].iloc[0]['runx1_pos']
                curdf["seqalign"] = curdf.apply(lambda x: align(x['Sequence'], s1-x['ets1_start']), axis=1)
                fig = plt.figure(figsize=(12,12))
                ds.plot_average(curdf, linemarks=[s1, s2], base_count=False, in_fig=fig, lblnames=["cooperative","independent"], pltlabel="orientation %s"%ori, label_title=True)
                pdf.savefig(fig)
                plt.close()
                with open("seqs.txt",'w') as f:
                    f.write("\n".join(curdf["seqalign"].tolist()))
