'''
Created on Oct 30, 2019

@author: vincentiusmartin
@editedby: Farica Zhuang
'''

import pandas as pd
import sys
import importlib

from chip2probe.util import util
from chip2probe.util import bio
from chip2probe.modeler.features import *

class CoopTrain:

    def __init__(self, trainingdata, corelen, sep="\t", flip_th=False, positive_cores=[]):
        """
        Cooperative Training class constructor

        Args:
            trainigdata : Path string or data frame with the required columns
            corelen : The length of the core binding site
            sep: Separator for the csv/tsv file if the trainingdata is string
            flip_th: Flip all TH orientation to HT (important whe we want to treat TH == HT)
            positive_cores: Required if flip_th is true, contain all the positive cores that define
            HT orientation.

         Returns:
            NA
        """
        if isinstance(trainingdata, pd.DataFrame):
            self.df = trainingdata.copy().reset_index(drop=True)
        elif isinstance(trainingdata, str):
            self.df = pd.read_csv(trainingpath, sep=sep)
        else:
            raise Exception("input must be string path or data frame")
        self.validate_cols()
        self.corelen = corelen
        if flip_th:
            if not positive_cores:
                raise Exception("Positive cores must be defined when flip_th is True")
            self.df = self.flip_th2ht(positive_cores)

    def validate_cols(self):
        """
        Check if we have all required columns.
        """
        required_cols = {"label"}
        if not required_cols.issubset(self.df.columns):
            raise Exception("Missing columns, required columns are: {}".format(required_cols))

    def flip_th2ht(self,positive_cores):
        # flip if orientation one
        ori = self.get_feature("orientation", {"positive_cores":positive_cores, "relative":True, "one_hot":False})
        records = self.df.to_dict(orient='records')
        for i in range(0,len(records)):
            if ori[i]["ori"] == "HT/TH":
                site_str = records[i]["sequence"][records[i]["site_str_pos"] - self.corelen//2:records[i]["site_str_pos"] + self.corelen//2]
                site_wk = records[i]["sequence"][records[i]["site_wk_pos"] - self.corelen//2:records[i]["site_wk_pos"] + self.corelen//2]
                if site_str not in positive_cores and site_wk not in positive_cores:
                    records[i]["sequence"] = bio.revcompstr(records[i]["sequence"])
                    # flip the position as well
                    str_pos = records[i]["site_str_pos"]
                    wk_pos = records[i]["site_wk_pos"]
                    records[i]["site_str_pos"] = len(records[i]["sequence"]) - wk_pos
                    records[i]["site_wk_pos"] = len(records[i]["sequence"]) - str_pos
        return pd.DataFrame(records)

    # ---------- GETTING FEATURE ----------
    def get_feature(self, feature, params):
        """
        Get a feature

        Args:
            feature: a string, the name of the feature (e.g. "orientation")
         Returns:
            feature_dict
        """

        module = importlib.import_module("chip2probe.modeler.features.%s" % feature)
        class_ = getattr(module, feature.capitalize())
        instance = class_(self.df, params)
        return instance.get_feature()

    def get_feature_all(self, feature_dict):
        """
        Get all feature based on feature dict

        Accepts a dictionary of feature name and parameters relative to the feature.
        All features are from the ``chip2probe.modeler.features`` package.

        Args:
            feature_dict: the following is the list of currently available feature:
                1. distance: {type:"numeric/categorical"}
                2. orientation: {"positive_cores:[]", relative:Bool, one_hot:Bool}
                3. affinity: {imads: sitespredict.imads instance}
         Returns:
            ldict: list of dictionary of features
        """
        ldict = []
        for class_name, params in feature_dict.items():
            ldict = util.merge_listdict(self.get_feature(class_name,params),ldict)
        #print(ldict)
        return ldict

    def get_numeric_label(self, label_map):
        ytrain = self.df['label'].map(label_map)
        return ytrain

    def get_training_df(self, feature_dict, label_map={'cooperative': 1, 'additive': 0}):
        """Get training df from a dictionary of features."""
        ldict = self.get_feature_all(feature_dict)
        train = pd.DataFrame(ldict)
        train['label'] = self.get_numeric_label(label_map).values
        return train
