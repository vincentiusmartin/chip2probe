import pickle
from chip2probe.modeler.cooptrain import CoopTrain
import pandas as pd

class ShapeModel:
    def __init__(self, paramdict):
        self.model = {}
        self.param = {}
        for key, val in paramdict.items():
            self.model[key] = pickle.load(open(val["path"], "rb"))
            self.param[key] = val["param"]

    def predict(self, df):
        # hardcoded for now, too tired
        ct = CoopTrain(df, corelen=4, flip_th=True, positive_cores=["GGAA","GGAT"])
        ori = ct.get_feature("orientation", {"positive_cores":["GGAA", "GGAT"]})
        df["orientation"] = pd.DataFrame(ori)["ori"]

        for key in self.model:
            curdf = df.loc[df["orientation"] == key]
            feature_dict = {
                    "distance":{"type":"numerical"},
                    "shape_in": {"seqin":4, "smode":"positional", "direction":"inout"}, # maximum seqin is 4
                    "shape_out": {"seqin":-4, "smode":"positional", "direction":"inout"}
                }
            print(key,self.param[key])
            train_df = pd.DataFrame(ct.get_feature_all(feature_dict))[self.param[key]]
            train = train_df.values.tolist()
            pred = self.model[key].predict(train)
            print(pred)
