import os
os.chdir("../../../..")
import pandas as pd
import pickle

import chip2probe.modeler.mutation as mut
import chip2probe.training_gen.traingen as tg
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.shapemodel import ShapeModel

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pbmescore import PBMEscore
from chip2probe.sitespredict.dnasequence import DNASequence

from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.kompas import Kompas

pd.set_option("display.max_columns",None)
if __name__ == "__main__":
    trainingpath = "output/array_design_files/Coop3Ets_validation/train.tsv"
    df = pd.read_csv(trainingpath,sep="\t")
    # select only genomic sequences
    df = df[~df['name'].str.contains("dist|weak")] # should already be all genomics

    # Load escore object
    escore = PBMEscore("input/site_models/escores/Ets1_8mers_11111111.txt")

    # Load imads object
    imads12_paths = ["input/site_models/imads_models/Ets1_w12_GGAA.model", "input/site_models/imads_models/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19) # 0.2128

    print("Number of input rows: %d"%df.shape[0])
    indf = pd.DataFrame(df[["sequence","label"]])# "id"
    indf["label"] = indf["label"].replace({"cooperative":1,"independent":0})

    mindist = imads12.sitewidth - 3
    ct = CoopTrain(indf["sequence"].values.tolist(), corelen=4, flip_th=True, imads=imads12, ignore_sites_err=True, seqcolname="sequence")
    om = ct.df.join(indf.set_index("sequence"), on="sequence", how="inner") # this already include the orientation
    seqs = []
    passval = 0
    for index, row in om.iterrows():
        sites, sites_specific = DNASequence(row["sequence"], imads12, escore, 0.4, 0).get_sites()
        if len(sites_specific) != 2: # or sites[1]["core_mid"] - sites[0]["core_mid"] < mindist: #or len(sites_specific) != 2
            #print(len(sites), sites[1]["core_mid"] - sites[0]["core_mid"], mindist)
           continue
        seqs.append(row["sequence"])
        passval += 1
    pd.DataFrame({'sequence':seqs}).to_csv("seqs.csv")
    print("Number of sequences passing the cutoff %d" % passval)

    # 1. Mutate based on affinity
    aff_m = mut.mutate_affinity(indf, imads12, escore, deep=6, idcol="name")
    print(aff_m["seqid"].unique())

    # 2. Mutate based on orientation
    ori_m = mut.mutate_orientation(indf, imads12,  escore, deep=6, idcol="name")

    # 3. Mutate based on distance
    dis_m = mut.mutate_dist(indf, imads12, escore, warning=False, deep=6, patch=False, idcol="name")
    print(dis_m[["seqid"]].drop_duplicates().shape[0])
    # dm.to_csv("custom_distance.csv", index=False, header=True)
    # dm = pd.read_csv("custom_distance.csv")
    # we need to add the orientation information for the distance

    mutdf = pd.concat([aff_m, dis_m, ori_m]) # ,
    mutdf.to_csv("custom.csv", index=False, header=True)
    mut_aff = mutdf.rename(columns={"sequence":"Sequence"})[["Sequence","site1_affinity","site2_affinity"]].drop_duplicates()

    pwm_ets = PWM("main_nar/input/sitemodels/ets1.txt", log=True)
    kompas_ets = Kompas("main_nar/input/sitemodels/Ets1_kmer_alignment.txt", core_start = 11, core_end = 15, core_center = 12)
    mutdf = pd.read_csv("custom.csv")[["seqid","sequence","muttype","comment","wtlabel"]].rename(columns={"sequence":"Sequence"})
    mutr = tg.gen_training(mutdf["Sequence"].tolist(), pwm_ets, kompas_ets)
    mutdf = mutdf.merge(mutr,on="Sequence")

    # not the most effective way since we calculate cooptrain twice...
    ct = CoopTrain(mutdf)
    feature_dict = {
        "distance":{"type":"numerical"},
        "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}},
        "shape_in":{"seqin":5, "poscols":['site_str_pos','site_wk_pos'], "smode":"relative"},
        "shape_out":{"seqin":-2, "poscols":['site_str_pos','site_wk_pos'], "smode":"relative"}
    }
    # cols = ['dist_numeric', 'ProT_outer_str_pos_2', 'MGW_inner_str_pos_3', 'ProT_outer_wk_pos_2', 'ProT_outer_wk_pos_1', 'Roll_inner_str_pos_2',  'ProT_inner_wk_pos_1', 'HelT_inner_str_pos_2', 'Roll_inner_wk_pos_2', 'MGW_inner_str_pos_2']
    train1 = ct.get_feature_all(feature_dict) #[cols]
    model1 = pickle.load(open("main_nar/output/Ets1Ets1/model/ets1_ets1_rfshapemodel.sav", "rb"))
    pred1 = model1.predict(train1)
    prob1 = model1.predict_proba(train1)
    mutdf["shape_pred"] = pred1
    mutdf["shape_prob"] = [prob1[i][1] for i in range(len(pred1))]

    feature_dict = {
        "distance":{"type":"numerical"},
        "affinity": {"colnames": ("site_str_score","site_wk_score")},
        "orientation": {"relative":True, "one_hot":True, "pos_cols": {"site_str_pos":"site_str_ori", "site_wk_pos":"site_wk_ori"}}
    }
    train = ct.get_feature_all(feature_dict)
    model = pickle.load(open("main_nar/output/Ets1Ets1/model/ets1_ets1_rfmodel.sav", "rb"))
    pred = model.predict(train)
    prob = model.predict_proba(train)
    mutdf["main_pred"] = pred
    mutdf["main_prob"] = [prob[i][1] for i in range(len(pred))] # use the probability of being cooperative
    seqnmap = df[["name","sequence"]].rename(columns={"sequence":"Sequence"})
    mutdf = mutdf.merge(seqnmap, on="Sequence", how="left")
    mutdf.merge(mut_aff, on="Sequence").to_csv("custom_withpred.csv", index=False, header=True)

    mutdf = pd.read_csv("custom_withpred.csv")

    # ------ Get statistics of the wild type sequences ------
    wtdf = mutdf.loc[mutdf["comment"] == "wt"][["seqid","Sequence","wtlabel","shape_pred","main_pred"]].drop_duplicates()
    main_true = wtdf.loc[wtdf["main_pred"] == wtdf["wtlabel"]]
    shape_true = wtdf.loc[wtdf["shape_pred"] == wtdf["wtlabel"]]
    intersect = main_true.merge(shape_true, on=["seqid","Sequence"])
    print("Number of wt sequence: %d" % wtdf.shape[0])
    print(("Number of correctly predicted wt training:\n" +
          "  main_pred : %d\n" +
          "  shape_pred: %d") % (main_true.shape[0], shape_true.shape[0]))
    print("Number of correctly predicted wt intersection main vs shape: %d" % intersect.shape[0])
