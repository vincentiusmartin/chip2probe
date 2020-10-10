import pandas as pd
import chip2probe.modeler.dnashape as ds
from chip2probe.util import bio

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
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1"
    df = pd.read_csv("%s/ch1_ch2/training_pwm.tsv" % basepath, sep="\t")

    fastadict = dict(zip(df["Name"], df["sequence"]))
    ds = ds.DNAShape(fastadict)

    curdf = pd.DataFrame(df[df['distance'] == 5])
    s1 = int(curdf['ets_pos'].mode())
    s2 = curdf[curdf['ets_pos'] == s1].iloc[0]['runx_pos']
    curdf["seqalign"] = curdf.apply(lambda x: align(x['sequence'], s1-x['ets_pos']), axis=1)

    ds.plot_average(curdf, linemarks=[s1, s2])
