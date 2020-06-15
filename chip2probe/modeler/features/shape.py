from chip2probe.modeler.features import basefeature
import chip2probe.util.bio as bio
import chip2probe.modeler.dnashape as ds
import string, random

from chip2probe.modeler.features.orientation import Orientation

import numpy as np

class Shape(basefeature.BaseFeature):
    def __init__(self, traindf, params):
        """
        DNA Shape feature prediction class

        Args:
            traindf: dataframe containing the "name", "sequence" column
            params:
                - c: trainingdata.dnashape.DNAShape object
                - self.seqin: if positive, get shape with direction to the inside;
                        if negative, get shape with direction to the outside
                - smode: "positional" or "strength" get direction anchored on the
                        position or strength of the site
                - direction= "inout" or "orientation" using in-out direction or based on
                        orientation
                - positive_cores = must be supplied if direction == "orientation"

         Returns:
            NA
        """
        default_args = {
            "seqin": 0,
            "smode": "positional", #site mode
            "direction": "inout",
            "positive_cores" : []
        }
        self.df = traindf
        self.set_attrs(params, default_args)
        if self.smode != "strength" and self.smode != "positional":
            raise TypeError("Smode can only be 'strength' or 'positional'")
        if self.direction != "inout" and self.direction != "orientation":
            raise TypeError("Direction can only be 'inout' or 'orientation'")
        if self.direction == "orientation" and ("positive_cores" not in params or not params["positive_cores"]):
            raise TypeError("Positive cores are needed when direction is 'orientation'")

        fastadict = dict(zip(self.df["name"], self.df["sequence"]))
        shapeobj = ds.DNAShape(fastadict)
        self.shapes = {k:getattr(shapeobj,k.lower()) for k in shapeobj.shapetypes}
        # make a dictionary of list instead of nested dictionary since we use this
        # as features
        namelist = self.df["name"].tolist()
        self.shapes = {k:[v[n] for n in namelist] for k, v in self.shapes.items()}
        if self.direction == "orientation":
            ori = Orientation(self.df, {"positive_cores":self.positive_cores}).get_feature()
            self.df["orientation"] = [o["ori"] for o in ori]

    def get_feature(self):
        """
        Get the shape features based on the features:
            - self.seqin:
            - smode: this determine whether we use the strength of binding sites
                as the anchor (i.e. weak, strong) of flank or the position (i.e.
                first site, second site)
            - direction:
        """
        rfeature = []
        for idx,row in self.df.iterrows():
            # first get which site is on the left (site1) and on the right (site2)
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
                if self.smode == "strength":
                    s1type, s2type = "str", "wk"
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
                if self.smode == "strength":
                    s1type, s2type = "wk", "str"
            # if site mode is positional, we use the position instead
            if self.smode == "positional":
                s1type, s2type = "s1", "s2"
            # get flanking shape based on direction
            if self.direction == "inout":
                rfeature.append(self.get_shape_inout(idx, site1, site2, type))
            else:
                rfeature.append(self.get_shape_orientation(idx, site1, site2, s1type, s2type , row["orientation"]))
        return rfeature

    def get_shape_inout(self, idx, site1, site2,  s1type, s2type):
        dfeature = {}
        for s in self.shapes:
            # orientation
            if self.seqin > 0: # inner
                flank1 = self.shapes[s][idx][site1:site1+self.seqin]
                flank2 = self.shapes[s][idx][site2-self.seqin:site2][::-1]
                type = "inner"
            else: # outer
                flank1 = self.shapes[s][idx][site1+self.seqin:site1][::-1]
                flank2 = self.shapes[s][idx][site2:site2-self.seqin]
                type = "outer"
            for i in range(abs(self.seqin)):
                dfeature["%s_%s_%s_pos_%d" % (s,type,s1type,i)] = flank1[i]
                dfeature["%s_%s_%s_pos_%d" % (s,type,s2type,i)] = flank2[i]
        return dfeature

    def get_shape_orientation(self, idx, site1, site2, s1type, s2type, orientation):
        dfeature = {}
        for s in self.shapes:
            # get the inner flanking region
            if orientation == 'HH':
                if self.seqin >= 0:
                    flank1 = self.shapes[s][idx][site1:site1 + self.seqin]
                    flank2 = self.shapes[s][idx][site2 - self.seqin:site2][::-1]
                    type = "head"
                if self.seqin < 0:
                    flank1 = self.shapes[s][idx][site1 + self.seqin:site1][::-1]
                    flank2 = self.shapes[s][idx][site2:site2 - self.seqin]
                    type = "tail"
            elif orientation == 'TT':
                if self.seqin < 0:
                    flank1 = self.shapes[s][idx][site1:site1 - self.seqin]
                    flank2 = self.shapes[s][idx][site2 + self.seqin:site2][::-1]
                    type = "tail"
                if self.seqin >= 0:
                    flank1 = self.shapes[s][idx][site1 - self.seqin:site1][::-1]
                    flank2 = self.shapes[s][idx][site2:site2 + self.seqin]
                    type = "head"
            elif orientation == 'HT/TH':
                if self.seqin >= 0:
                    flank1 = self.shapes[s][idx][site1:site1 + self.seqin]
                    flank2 = self.shapes[s][idx][site2:site2 + self.seqin]
                    type = "head"
                if self.seqin < 0:
                    flank1 = self.shapes[s][idx][site1 + self.seqin:site1][::-1]
                    flank2 = self.shapes[s][idx][site2 + self.seqin:site2][::-1]
                    type = "tail"
            for i in range(abs(self.seqin)):
                dfeature["%s_%s_%s_pos_%d" % (s,type,s1type,i)] = flank1[i]
                dfeature["%s_%s_%s_pos_%d" % (s,type,s2type,i)] = flank2[i]
        return dfeature
