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
    pd.set_option("display.max_rows",None)

    # basepath = "output/Ets1Runx1"
    # df = pd.read_csv("%s/training/train_ets1_runx1.tsv"%basepath, sep="\t")
    # poscol1 = "ets1_pos"
    # poscol2 = "runx1_pos"

    basepath = "output/Ets1Ets1"
    df = pd.read_csv("%s/training/train_ets1_ets1.tsv"%basepath, sep="\t")
    poscol1 = "site_str_pos"
    poscol2 = "site_wk_pos"

    df["s1pos"] = df.apply(lambda x: x[poscol1] if x[poscol1] < x[poscol2] else x[poscol2], axis=1)
    df["s2pos"] = df.apply(lambda x: x[poscol1] if x[poscol1] > x[poscol2] else x[poscol2], axis=1)


    dist = sorted(df["distance"].unique().tolist())[1:]
    oris = df["orientation"].unique().tolist()

    fastadict = dict(zip(df["Name"], df["Sequence"]))
    ds = ds.DNAShape(fastadict)
    n = 0
    labels = ["independent", "cooperative"]

    print(df.groupby(["distance","orientation"]).count()[["Name"]])

    with PdfPages("%s/model/shapepos/shape_plots.pdf" % basepath) as pdf:
        for d in dist:
            curdf_dist = pd.DataFrame(df[df['distance'] == d])
            for ori in oris: # oris
                curdf = pd.DataFrame(curdf_dist[curdf_dist['orientation'] == ori])
                s1 = int(curdf['s1pos'].mode().iloc[0])
                s2 = curdf[curdf['s1pos'] == s1].iloc[0]['s2pos']
                curdf["seqalign"] = curdf.apply(lambda x: align(x['Sequence'], s1-x['s1pos']), axis=1)
                fig = plt.figure(figsize=(12,12))
                ds.plot_average(curdf, linemarks=[s1, s2], base_count=False, in_fig=fig, lblnames=["cooperative","independent"], pltlabel="orientation %s"%ori, label_title=True)
                pdf.savefig(fig)
                plt.close()
                oristr = ori.replace("/","")
                for l in labels:
                    curdflbl = curdf[curdf["label"] == l]
                    with open("%s/model/shapepos/seqs_dist%d_%s_%s.txt"%(basepath,d,oristr,l),'w') as f:
                        f.write("\n".join(curdflbl["seqalign"].tolist()))
