'''
Created on Jul 16, 2019

@author: vincentiusmartin
'''

import sys

import abc
import pandas as pd
import util.bio as bio

class BaseModel(abc.ABC):
    '''
    classdocs
    '''

    @abc.abstractclassmethod
    def __init__(self):
        pass

    def pred_input_todict(self, sequence_input, sequence_colname="sequence", key_colname=""):
        if isinstance(sequence_input, pd.DataFrame):
            return bio.get_seqdict(sequence_input, sequence_colname=sequence_colname, keycolname=key_colname)
        elif isinstance(sequence_input, dict):
            return sequence_input
        else:
            raise Exception("input must be data frame or dictionary of sequences")

    @abc.abstractclassmethod
    def predict_sequences(self, sequence_df):
        """
        Should return basepred
        """
        pass

    @abc.abstractclassmethod
    def make_plot_data(self, sequence_df):
        pass
