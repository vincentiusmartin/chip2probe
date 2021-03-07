import os
os.chdir("../../../..")
import pandas as pd
import pickle

import chip2probe.modeler.mutation as mut
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.shapemodel import ShapeModel

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.pbmescore import PBMEscore
from chip2probe.sitespredict.dnasequence import DNASequence

if __name__ == "__main__":
    trainingpath = "output/homotypic/training/training.csv"
    df = pd.read_csv(trainingpath)
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
    indf["label"] = indf["label"].replace({"cooperative":1,"additive":0})

    mindist = imads12.sitewidth - 3
    ct = CoopTrain(indf["sequence"].values.tolist(), corelen=4, flip_th=True, imads=imads12, ignore_sites_err=True)
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

    """
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
    mutdf = pd.read_csv("custom.csv")

    # not the most effective way since we calculate cooptrain twice...
    ct = CoopTrain(mutdf["sequence"].values.tolist(), corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"], imads=imads12)
    feature_dict = {
        "distance":{"type":"numerical"},
         "shape_in": {"seqin":4, "smode":"strength", "direction":"inout"}, # maximum seqin is 4
         "shape_out": {"seqin":-4, "smode":"strength", "direction":"inout"}
    }
    cols = ['dist_numeric', 'ProT_outer_str_pos_2', 'MGW_inner_str_pos_3', 'ProT_outer_wk_pos_2', 'ProT_outer_wk_pos_1', 'Roll_inner_str_pos_2',  'ProT_inner_wk_pos_1', 'HelT_inner_str_pos_2', 'Roll_inner_wk_pos_2', 'MGW_inner_str_pos_2']
    train_df1 = pd.DataFrame(ct.get_feature_all(feature_dict))[cols]
    train1 = train_df1.values.tolist()
    model1 = pickle.load(open("input/modeler/coopmodel/dist_shapeio.sav", "rb"))
    pred1 = model1.predict(train1)
    prob1 = model1.predict_proba(train1)
    mutdf["shape_pred"] = pred1
    mutdf["shape_prob"] = [prob1[i][1] for i in range(len(pred1))]

    feature_dict = {
        "distance":{"type":"numerical"},
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True},
        "affinity": {"imads":imads12}
    }
    train = ct.get_feature_all(feature_dict, rettype="list")
    model = pickle.load(open("input/modeler/coopmodel/dist_ori_12merimads.sav", "rb"))
    pred = model.predict(train)
    prob = model.predict_proba(train)
    mutdf["main_pred"] = pred
    mutdf["main_prob"] = [prob[i][1] for i in range(len(pred))] # use the probability of being cooperative
    mutdf.to_csv("custom_withpred.csv", index=False, header=True)

    mutdf = pd.read_csv("custom_withpred.csv")

    # ------ Get statistics of the wild type sequences ------
    wtdf = mutdf.loc[mutdf["comment"] == "wt"][["seqid","sequence","wtlabel","shape_pred","main_pred"]].drop_duplicates()
    main_true = wtdf.loc[wtdf["main_pred"] == wtdf["wtlabel"]]
    shape_true = wtdf.loc[wtdf["shape_pred"] == wtdf["wtlabel"]]
    intersect = main_true.merge(shape_true, on=["seqid","sequence"])
    print("Number of wt sequence: %d" % wtdf.shape[0])
    print(("Number of correctly predicted wt training:\n" +
          "  main_pred : %d\n" +
          "  shape_pred: %d") % (main_true.shape[0], shape_true.shape[0]))
    print("Number of correctly predicted wt intersection main vs shape: %d" % intersect.shape[0])
    """
