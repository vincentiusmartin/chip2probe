'''
Created on Aug 21, 2019

@author: vincentiusmartin
'''

class CustomSeqMaker(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
        pass
    
    def eliminate_site(self, sequence, pbmescore, escore_threshold = 0.35):
        """
        This assumes that input is a sequence of a site
        """
        maxescore = 1 # just to initialize with a large value
        prevmax = -1 # just to avoid infinite loop
        site_mutpos = []
        full_sequence = str(sequence)
        
        flag = True
        while flag and prevmax != maxescore:    
            # INITIALIZE max escore calculation
            # if we don't need non-core-intersecting: 
            #epreds = pbmescore.predict_sequence(full_sequence)
            #maxepreds = max(epreds.predictions, key=lambda x:x['score'])
            prevmax = float(maxescore)
            maxepreds = self.get_non_core_intersecting_maxescore(full_sequence, site_index, pbmescore)
            maxescore = maxepreds['score'] 
            if maxescore < escore_threshold: # our conidition is met
                flag = False
            else:
                if maxescore == float("inf"): # no e-score that can be chosen
                    print("No e-score site can be mutated for sequence %s" % full_sequence)
                    # since there is no site to be mutated, then just use empty list
                    return full_sequence, []
                seq_tomutate = maxepreds["escore_seq"]
                midpos = len(seq_tomutate) // 2
                mutated_escore_seq = self.mutate_escore_seq_at_pos(seq_tomutate, midpos, pbmescore)
                
                # mutate the sequence
                mut_start = self.bsites[site_index].site_start + maxepreds["start_idx"]
                mut_end = mut_start + len(mutated_escore_seq)
                full_sequence = full_sequence[:mut_start] + mutated_escore_seq + full_sequence[mut_end:]
                site_mutpos.append(mut_start + midpos)
        return full_sequence, site_mutpos