import pandas as pd

import chip2probe.training_gen.arranalysis as arr

if __name__=="__main__":
    pd.set_option("display.max_columns",None)
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/df_labeled/ch1_vs_ch2"
    df1, df2 = pd.read_csv("%s/probes_labeled_er.csv"%basepath), pd.read_csv("%s/probes_labeled_re.csv"%basepath)
    arr.plot_ori_inconsistency(df1, df2, cols=["ch1","ch2"],orinames=["er","re"],log=True)
