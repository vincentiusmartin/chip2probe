import os
os.chdir("../..")
import pandas as pd
import pickle

import chip2probe.modeler.mutation as mut
from chip2probe.probe_generator.probefilter.sitespredict.imads import iMADS
from chip2probe.probe_generator.probefilter.sitespredict.imadsmodel import iMADSModel
from chip2probe.modeler.cooptrain import CoopTrain
from chip2probe.modeler.shapemodel import ShapeModel

if __name__ == "__main__":
    trainingpath = "input/modeler/training_data/training_p01_adjusted.tsv"
    df = pd.read_csv(trainingpath, sep="\t")
    # select only genomic sequences
    df = df[~df['name'].str.contains("dist|weak")]

    imads12_paths = ["input/modeler/imads_model/Ets1_w12_GGAA.model", "input/modeler/imads_model/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19) # 0.2128

    print("Number of input rows: %d"%df.shape[0])
    indf = df[["id","sequence","label"]]
    indf["label"] = indf["label"].replace({"cooperative":1,"additive":0})

    # 1. Mutate based on affinity
    aff_m = mut.mutate_affinity(indf, imads12, deep=3)

    # 2. Mutate based on orientation
    ori_m = mut.mutate_orientation(indf, imads12, deep=3)

    # 3. Mutate based on distance
    dis_m = mut.mutate_dist(indf, imads12, warning=False, deep=3)
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
    train = ct.get_feature_all(feature_dict, aslist=True)
    model = pickle.load(open("input/modeler/coopmodel/dist_ori_12merimads.sav", "rb"))
    pred = model.predict(train)
    prob = model.predict_proba(train)
    mutdf["main_pred"] = pred
    mutdf["main_prob"] = [prob[i][1] for i in range(len(pred))] # use the probability of being cooperative
    mutdf.to_csv("custom_withpred.csv", index=False, header=True)

    mutdf = pd.read_csv("custom_withpred.csv")
    wtdf = mutdf.loc[mutdf["comment"] == "wt"][["seqid","sequence","wtlabel","shape_pred","main_pred"]].drop_duplicates()
    main_true = wtdf.loc[wtdf["main_pred"] == wtdf["wtlabel"]]
    shape_true = wtdf.loc[wtdf["shape_pred"] == wtdf["wtlabel"]]
    intersect = main_true.merge(shape_true, on=["seqid","sequence"])
    print(main_true.shape[0], shape_true.shape[0])
    print("Number of wt sequence: %d" % wtdf.shape[0])
    print(("Number of correctly predicted wt training:\n" +
          "  main_pred : %d\n" +
          "  shape_pred: %d") % (main_true.shape[0], shape_true.shape[0]))
    print("Number of correctly predicted wt intersection main vs shape: %d" % intersect.shape[0])
