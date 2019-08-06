'''
Created on Jul 23, 2019

@author: vincentiusmartin
'''
from sitesfinder.sitesfinder import SitesFinder

# UNDER CONSTRUCTION
class PWM(SitesFinder):
    '''
    classdocs
    '''


    def __init__(self, sequence_df, pwm_path, ):
        super(self.__class__, self).__init__(sequence_df)
        
    def read_dna_pwm(self, pwm_path):
        """
        Read DNA (i.e. A,C,G,T) PWM
        """
        pass
        
    def predict_sequences(self):
        pass
    
    def plot(self, ax):
        pass
        