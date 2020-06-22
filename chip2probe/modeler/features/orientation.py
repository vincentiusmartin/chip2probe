import pandas as pd

from chip2probe.util import bio
from chip2probe.modeler.features import basefeature

class Orientation(basefeature.BaseFeature):
    def __init__(self, traindf, params):
        """
        Orientation feature class

        Args:
            traindf: dataframe containing the sequence column
            params:
                - positive_cores: list of string of positive cores
                - relative: use site orientation relative to each other (HH, HT/TH, TT)
                - one_hot: use one_hot representation

         Returns:
            NA
        """

        # Make positive cores as mandatory
        if "positive_cores" not in params or not params["positive_cores"]:
            raise Exception("Positive cores are needed")

        default_args = {
            "positive_cores":[],
            "relative": True,
            "one_hot": False
        }
        self.df = traindf
        self.set_attrs(params, default_args)
        self.motiflen = len(self.positive_cores[0])

    def get_feature(self):
        """

        """
        negative_cores = [bio.revcompstr(p) for p in self.positive_cores]
        rfeature = []
        for idx,row in self.df.iterrows():
            # get the binding sites
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
            seq = row["sequence"]
            p1 = seq[site1 - self.motiflen//2:site1 + self.motiflen//2]
            p2 = seq[site2 - self.motiflen//2:site2 + self.motiflen//2]
            if p1 in self.positive_cores:
                s1 = 1
            elif p1 in negative_cores:
                s1 = -1
            else:
                s1 = 0
                print("couldn't find the first site %s in %s in the core list" % (p1,seq))

            if p2 in self.positive_cores:
                s2 = 1
            elif p2 in negative_cores:
                s2 = -1
            else:
                s2 = 0
                print("couldn't find the second site %s in %s in the core list" % (p2,seq))
            if self.relative:
                if s1 == 1 and s2 == 1:
                    ori = 'HT/TH'
                elif s1 == -1 and s2 == -1:
                    ori = 'HT/TH'
                elif s1 == 1 and s2 == -1:
                    ori = 'HH'
                elif s1 == -1 and s2 == 1:
                    ori = 'TT'
                else:
                    ori = '-1'
                rfeature.append({"ori":ori})
            else:
                rfeature.append({"ori1":s1, "ori2":s2})
        if self.relative and self.one_hot:
            dum_df = pd.DataFrame(rfeature)
            notfound = {"HH","HT/TH","TT"} - set(dum_df["ori"].unique())
            dum_rec = pd.get_dummies(dum_df).to_dict('records')
            for nf in notfound:
                print("notfound orientation in the dataset: ",nf)
                for i in range(len(dum_rec)):
                    dum_rec[i]["ori_%s" % nf] = 0
            return dum_rec
        else:
            return rfeature
