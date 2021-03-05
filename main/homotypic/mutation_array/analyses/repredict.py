import pandas as pd
import os
import pickle
from chip2probe.modeler.cooptrain import CoopTrain

from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel

os.chdir("../../../..")

if __name__ == "__main__":
    imads12_paths = ["input/site_models/imads_models/Ets1_w12_GGAA.model", "input/site_models/imads_models/Ets1_w12_GGAT.model"]
    imads12_cores = ["GGAA", "GGAT"]
    imads12_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads12_paths, imads12_cores)]
    imads12 = iMADS(imads12_models, 0.19) # 0.2128

    selected = pd.read_csv("output/array_design_files/Coop2Ets_validation/custom_probes_selected.csv")
    selectedlist = selected["sequence"].values.tolist()
    model = pickle.load(open("input/modeler/coopmodel/dist_ori_12merimads.sav", "rb"))
    ct = CoopTrain(selectedlist, corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"], imads=imads12)
    feature_dict = {
        "distance":{"type":"numerical"},
        "orientation": {"positive_cores":["GGAA","GGAT"], "one_hot":True},
        "affinity": {"imads":imads12}
    }
    train = pd.DataFrame(ct.get_feature_all(feature_dict)).values.tolist()
    # pickle.dump(train, open("train.pickle", 'wb'))
    # train = pickle.load(open("train.pickle", "rb"))
    pred = model.predict(train)
    prob = model.predict_proba(train)
    selected["main_pred"] = pred
    selected["main_prob"] = [prob[i][1] for i in range(len(prob))]
    selected.to_csv("custom_probes_selected_repred.csv", index=False)
