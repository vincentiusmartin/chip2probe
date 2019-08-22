'''
Created on Jul 22, 2019

@author: vincentiusmartin
'''
from sitesfinder.sitesfinder import SitesFinder
import util.bio as bio 
from sitesfinder.prediction.basepred import BasePrediction

class PBMEscore(SitesFinder):
    '''
    classdocs
    '''


    def __init__(self, escore_short_path, escore_map_path, kmer=8):
        '''
        Constructor
        '''
        self.escore = self.read_escore(escore_short_path, escore_map_path)
        self.kmer = kmer
        
    
    def read_escore(self, escore_short_path, escore_map_path):
        with open(escore_short_path) as f:
            eshort = [float(line) for line in f]
        with open(escore_map_path) as f:
            next(f)
            emap = [int(line.split(",")[1])-1 for line in f]
        elong = [eshort[idx] for idx in emap]
        return elong
            
    def predict_sequence(self, sequence):
        """
        input: sequence string
        """
        prediction = []
        for i in range(0,len(sequence)-self.kmer+1):
            score = self.escore[bio.seqtoi(sequence[i:i+self.kmer])]
            # list index is the position from the first e-scoreo (i.e. i-th position)
            prediction.append({"position":i+(self.kmer+1)//2,"escore_seq":sequence[i:i+self.kmer],"score":score,"start_idx":i})
        return BasePrediction(sequence, prediction)
    
    # TODO: modify sequence to use this function instead
    def get_escores_specific(self, sequence, escore_cutoff = 0.4, escore_gap = 0):
        escores = self.predict_sequence(sequence).predictions
        signifcount = 0
        startidx = -1
        gapcount = 0
        escore_signifsites = []
        for i in range(0, len(escores)):
            escoresite = escores[i]
            if escoresite["score"] > escore_cutoff :
                if signifcount == 0:
                    startidx = i
                signifcount += 1
                gapcount = 0 
            # we can ignore else if here since we need i == len(esores)-1
            if escoresite["score"] <= escore_cutoff and i != len(escores)-1 and gapcount < escore_gap:
                # check if the sequence is still within 
                gapcount += 1
            elif escoresite["score"] <= escore_cutoff or i == len(escores)-1: 
                if signifcount > 0:
                    # if we have found sufficient e-scores above the cutoff then get the binding sites
                    if signifcount >= 2:
                        # startpos: the start of binding
                        escore_bind = {"startpos":escores[startidx]['position'],  "escorelength":signifcount + gapcount, 
                                "escore_startidx":escores[startidx]['start_idx']}
                        escore_signifsites.append(escore_bind)
                    startidx = -1
                    signifcount = 0  
        return escore_signifsites
        
    
    def predict_sequences(self, sequence_df, sequence_colname = "sequence"):
        """
        input: sequence data frame or sequence dictionary
        """
        seqdict = self.pred_input_todict(sequence_df, sequence_colname=sequence_colname)
        predictions = {}
        for key,sequence in seqdict.items():
            predictions[key] = self.predict_sequence(sequence)
        return predictions
            
    def plot(self, predictions_dict, scale = 1, escore_cutoff = 0.4, additional_functions = {}):
        func_dict = {}
        for key in predictions_dict:
            sequence = predictions_dict[key].sequence
            escores = predictions_dict[key].predictions
            func_pred =[]
            y_escore = [x["score"] * scale for x in escores]
            x_escore = [x["position"] for x in escores]
            func_pred.append({"func":"plot","args":[x_escore, y_escore],"kwargs":{"color":"orange", "linewidth" : 2.5}})
            func_pred.append({"func":"axhline", "args":[escore_cutoff * scale], "kwargs":{"color":"darkorange", "linestyle" : "dashed", "linewidth":1}})
            if key in additional_functions and additional_functions[key]:
                func_pred.extend(additional_functions[key])
            func_dict[key] = {"sequence":sequence,
                                  "plt":func_pred}
        return func_dict
    