from chip2probe.modeler.cooptrain import CoopTrain
import chip2probe.training_gen.traingen as tg
import pandas as pd

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.kompas import Kompas

import chip2probe.modeler.plotlib as pl

def get_sites_pos(df, kompas, pwm, seqcol="Sequence"):
    """
    Get site position for each sequence

    Args:
        df: input data frame

    """
    seqlist = df[seqcol].unique().tolist()
    poslist = []
    misscount = 0
    for seq in seqlist:
        x = kompas.predict_sequence(seq)
        if len(x) != 2:
            continue
        poslist.append({seqcol:seq, "site_str_pos":x[0]['core_start']+2, 'site_str_start':x[0]['core_start'], 'site_wk_pos':x[1]['core_start']+2, 'site_wk_start':x[1]['core_start']})
    posdf = pd.DataFrame(poslist)
    posdf['site_str_score'], posdf['site_str_ori'], posdf['site_str_core'] =  tg.pwm_score(posdf, pwm, "site_str_start", 4, 3, seqcol=seqcol)
    posdf['site_wk_score'], posdf['site_wk_ori'],  posdf['site_wk_core'] =  tg.pwm_score(posdf, pwm, "site_wk_start", 4, 3, seqcol=seqcol)

    # currently site str is the first site
    posdf['site_str_ori'].replace({0:"-",1:"+"}, inplace=True)
    posdf['site_wk_ori'].replace({0:"-",1:"+"}, inplace=True)
    posdf["orientation"] = posdf.apply(lambda x: "%s/%s" % (x["site_str_ori"],x["site_wk_ori"]), axis=1)
    posdf["distance"] = posdf["site_wk_pos"] - posdf["site_str_pos"]

    flip_target = []
    for i,r in posdf.iterrows():
        if r["site_str_score"] < r["site_wk_score"]:
            flip_target.append(i)

    posdf.loc[flip_target,['site_str_score','site_wk_score']] = posdf.loc[flip_target,['site_wk_score','site_str_score']].values
    posdf.loc[flip_target,['site_str_ori','site_wk_ori']] = posdf.loc[flip_target,['site_wk_ori','site_str_ori']].values
    posdf.loc[flip_target,['site_str_core','site_wk_core']] = posdf.loc[flip_target,['site_wk_core','site_str_core']].values
    posdf = posdf[[seqcol,"site_str_pos","site_str_score","site_wk_pos","site_wk_score" ,"distance","orientation"]] #
    return df.merge(posdf,on=seqcol) # "site_str_ori","site_str_core", "site_wk_ori","site_wk_core"

if __name__ == "__main__":
    pd.set_option("display.max_columns",None)

    # using pwm
    pwm_ets = PWM("input/sitemodels/ets1.txt", log=True)
    kompas_ets = Kompas("input/sitemodels/Ets1_kmer_alignment.txt",
                    core_start = 11, core_end = 15, core_center = 12)
    df = pd.read_csv("output/Ets1Ets1/label_pr/ets_ets_seqlabeled.csv").drop_duplicates()
    df = df[(df["label"] == "cooperative") | (df["label"] == "independent")]

    train = get_sites_pos(df, kompas_ets, pwm_ets)
    train = train[(train["site_wk_score"] != - 999) & (train["site_str_score"] != - 999)]
    train["orientation"].replace({"-/+":"+/-"}, inplace=True)
    print(train["label"].value_counts())
    train.to_csv("train_ets_ets.csv", index=False)

    train.rename(columns={'site_str_score': 'Binding strength of the stronger site', 'site_wk_score': 'Binding strength of the weaker site'}, inplace=True)
    pl.plot_stacked_categories(train, "distance", path="distance_bar.png", title="Distance distribution", ratio=True, figsize=(17,4))
    pl.plot_stacked_categories(train, "orientation", path="ori_bar.png", title="Relative sites orientation\ndistribution", ratio=True, figsize=(9,5))
    pl.plot_box_categories(train, incols=["Binding strength of the stronger site", "Binding strength of the weaker site"], alternative="smaller")
