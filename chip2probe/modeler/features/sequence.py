from chip2probe.modeler.features import basefeature
import chip2probe.util.bio as bio
import itertools

import numpy as np

class Sequence(basefeature.BaseFeature):
    def __init__(self, traindf, params):
        """
        DNA sequence feature prediction class

        Args:
            traindf: dataframe containing the "name", "sequence" column
            params:


         Returns:
            NA
        """
        default_args = {
            "seqin": 0,
            "poscols": [],
            "namecol":"name"
        }
        self.df = traindf
        self.set_attrs(params, default_args)

    def get_feature(self):
        rfeature = []
        for idx, row in self.df.iterrows():
            site1, site2 = row[self.poscols[0]], row[self.poscols[1]]
            flank1 = row["sequence"][site1:site1 + self.seqin] if self.seqin > 0 else row["sequence"][site1 + self.seqin:site1][::-1]
            flank2 = row["sequence"][site2 - self.seqin:site2][::-1] if self.seqin > 0 else row["sequence"][site2:site2 - self.seqin]
            label = "flank_in" if self.seqin > 0 else "flank_out"
            d1 = self.extract_positional(flank1, label=label)
            d2 = self.extract_positional(flank2, label=label)
            rfeature.append({**d1, **d2})
        return rfeature


    def extract_positional(self,seq, maxk=2, label="seq", minseqlen=-float("inf")):
        '''
        orientation: if right, then start from 0 to the right, else start from
        len(seq)-1 to the left
        minseqlen: will fill -1 if not enough bases
        '''
        iterseq = str(seq)
        nucleotides = ['A','C','G','T']
        features = {}
        for k in range(1,maxk+1):
            perm = ["".join(p) for p in itertools.product(nucleotides, repeat=k)]
            i = 0
            while i < len(iterseq)+1-k:
                for kmer in perm:
                    seqcmp = iterseq[i:i+k]
                    if seqcmp == kmer:
                        features["%s_pos%d_%s" % (label,i,kmer)] = 1
                    else:
                        features["%s_pos%d_%s" % (label,i,kmer)] = 0
                i += 1
            # append the rest with -1
            if minseqlen > 0:
                while i < minseqlen + 1 - k:
                    for kmer in perm:
                        features["%s_pos%d_%s" % (label,i,kmer)] = -1
                    i += 1
        return features
