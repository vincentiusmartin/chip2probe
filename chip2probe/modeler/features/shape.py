from chip2probe.modeler.features import basefeature
import chip2probe.util.bio as bio
import string, random

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

class Shape(basefeature.BaseFeature):
    def __init__(self, traindf, params):
        """
        DNA Shape feature prediction class

        Args:
            traindf: dataframe containing the "name", "sequence" column
            params:
                - c: trainingdata.dnashape.DNAShape object
                - seqin: if positive, get shape with direction to the inside;
                        if negative, get shape with direction to the outside
                - smode: "positional" or "strength" get direction anchored on the
                        position or strength of the site
                - direction= "inout" or "ori" using in-out direction or based on
                        orientation

         Returns:
            NA
        """
        default_args = {
            "seqin": 0,
            "smode": "positional",
            "direction": "inout"
        }
        self.df = traindf
        self.set_attrs(params, default_args)

        fastadict = dict(zip(self.df["name"], self.df["sequence"]))
        tmpfasta = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
        dnashape_r = importr('DNAshapeR')
        ds = dnashape_r.getShape("test.fasta")
        print(ds)
        #bio.makefasta(fastadict,"%s.fasta" % tmpfasta)
        import sys
        sys.exit(0)

    def get_feature(self):
        self.get_feature_flank_shapes()
        return None

    def get_feature_flank_shapes_orientation(self, dnashape, seqin, site_mode="strength"):
        if site_mode!="strength" and site_mode!="positional":
            raise Exception("Site mode can only be 'strength' or 'positional'")
        shapes = {"prot":dnashape.prot, "mgw":dnashape.mgw, "roll":dnashape.roll, "helt":dnashape.helt}
        rfeature = []
        for idx,row in self.df.iterrows():
            orientation = row["orientation"]
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
                if site_mode == "strength":
                    s1type, s2type = "str", "wk"
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
                if site_mode == "strength":
                    s1type, s2type = "wk", "str"
            if site_mode == "positional":
                s1type, s2type = "s1", "s2"
            dfeature = {}
            for s in shapes:
                # get the inner flanking region
                if orientation == 'HH':
                    if seqin >= 0:
                        flank1 = shapes[s][str(idx + 1)][site1:site1 + seqin]
                        flank2 = shapes[s][str(idx + 1)][site2 - seqin:site2][::-1]
                        type = "head"
                    if seqin < 0:
                        flank1 = shapes[s][str(idx + 1)][site1 + seqin:site1][::-1]
                        flank2 = shapes[s][str(idx + 1)][site2:site2 - seqin]
                        type = "tail"
                elif orientation == 'TT':
                    if seqin < 0:
                        flank1 = shapes[s][str(idx + 1)][site1:site1 - seqin]
                        flank2 = shapes[s][str(idx + 1)][site2 + seqin:site2][::-1]
                        type = "tail"
                    if seqin >= 0:
                        flank1 = shapes[s][str(idx + 1)][site1 - seqin:site1][::-1]
                        flank2 = shapes[s][str(idx + 1)][site2:site2 + seqin]
                        type = "head"
                elif orientation == 'HT/TH':
                    if seqin >= 0:
                        flank1 = shapes[s][str(idx + 1)][site1:site1 + seqin]
                        flank2 = shapes[s][str(idx + 1)][site2:site2 + seqin]
                        type = "head"
                    if seqin < 0:
                        flank1 = shapes[s][str(idx + 1)][site1 + seqin:site1][::-1]
                        flank2 = shapes[s][str(idx + 1)][site2 + seqin:site2][::-1]
                        type = "tail"
                for i in range(abs(seqin)):
                    dfeature["%s_%s_%s_pos_%d" % (s,type,s1type,i)] = flank1[i]
                    dfeature["%s_%s_%s_pos_%d" % (s,type,s2type,i)] = flank2[i]
            rfeature.append(dfeature)
        return rfeature

    def get_feature_flank_shapes(self):
        if self.smode != "strength" and self.smode != "positional":
            raise ValueError("Site mode can only be 'strength' or 'positional'")
        shapes = {"prot":self.dnashape.prot, "mgw":self.dnashape.mgw, "roll":self.dnashape.roll, "helt":self.dnashape.helt}
        rfeature = []
        for idx,row in self.df.iterrows():
            if row["site_wk_pos"] > row["site_str_pos"]:
                site1, site2 = row["site_str_pos"], row["site_wk_pos"]
                if site_mode == "strength":
                    s1type, s2type = "str", "wk"
            else:
                site1, site2 = row["site_wk_pos"], row["site_str_pos"]
                if site_mode == "strength":
                    s1type, s2type = "wk", "str"
            if site_mode == "positional":
                s1type, s2type = "s1", "s2"
            dfeature = {}
            for s in shapes:
                if seqin > 0: # inner
                    flank1 = shapes[s][str(idx + 1)][site1:site1+seqin]
                    flank2 = shapes[s][str(idx + 1)][site2-seqin:site2][::-1]
                    type = "inner"
                else: # outer
                    flank1 = shapes[s][str(idx + 1)][site1+seqin:site1][::-1]
                    flank2 = shapes[s][str(idx + 1)][site2:site2-seqin]
                    type = "outer"
                for i in range(abs(seqin)):
                    dfeature["%s_%s_%s_pos_%d" % (s,type,s1type,i)] = flank1[i]
                    dfeature["%s_%s_%s_pos_%d" % (s,type,s2type,i)] = flank2[i]
            rfeature.append(dfeature)
        return rfeature
