
from chip2probe.modeler.features import basefeature

class Affinity(basefeature.BaseFeature):
    def __init__(self, traindf, params):
        """
        Affinity prediction feature class

        Args:
            traindf: dataframe containing the site_wk_score and site_str_score columns
                     if imadsmodel is None
            params:
                - imadsmodel: imadsmodel,
                - modelwidth: width

         Returns:
            NA
        """
        default_args = {
            "imads" : None
        }
        self.df = traindf
        self.set_attrs(params, default_args)


    def get_feature(self):
        """Get a dictionary of binding site preference scores."""
        rfeature = []
        # if imadsmodel is not provided, we assume the weak and strong site position is provided
        # TODO: make general column name
        for idx, row in self.df.iterrows():
            if not self.imads:
                f = {"site_wk_score": row["site_wk_score"],
                     "site_str_score": row["site_str_score"]}
            else: # imads model is provided
                pr = self.imads.predict_sequence(row["sequence"])
                if len(pr) != 2:
                    #print("Found sites more than 2 sites, using site position in the table as reference")
                    newpr = []
                    # if we suddenly have > 2 predictions, try using position in the training
                    for p in pr:
                        if p["core_mid"] == row["site_str_pos"] or p["core_mid"] == row["site_wk_pos"]:
                            newpr.append(p)
                    if len(newpr) != 2:
                        print("Error on row:\n", row, "with prediction\n", pr)
                        raise Exception("Sequence does not have 2 binding sites")
                    pr = newpr
                f = {"site_wk_score": pr[0]["score"] if pr[0]["score"] < pr[1]["score"] else pr[1]["score"],
                     "site_str_score": pr[0]["score"] if pr[0]["score"] > pr[1]["score"] else pr[1]["score"]}
            rfeature.append(f)
        return rfeature
