import pandas as pd

# all permutations are already reverse-deleted
# all sequences are represented in binary


def revcompstr(seq):
    rev = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([rev[base] for base in reversed(seq)])


def get_seqdict(sequence_tbl, sequence_colname = "sequence", keycolname="" , keyname="sequence", 
                ignore_missing_colname=False):
    """
    Generates dictionary representation from sequence table or list
    :param sequence_tbl: 
    :param sequence_colname: 
    :param ignore_missing_colname: 
    :return: 
    """
    if not isinstance(sequence_tbl, pd.DataFrame) and not isinstance(sequence_tbl, list):  # and not isinstance(sequence_tbl, dict):
        raise Exception('sequence_tbl must be either pandas data frame or a list or a dictionary')
    if type(sequence_tbl) == pd.DataFrame:
        if not sequence_colname in sequence_tbl:
            if ignore_missing_colname:
                seqlist = [""] * sequence_tbl.shape[0]
            else:
                raise Exception('Could not find column %s in the data frame. This can be ignored by turning on \'ignore_missing_colname\'.' % sequence_colname)
        else:
            if keycolname: # using key from one of the column
                return pd.Series(sequence_tbl[sequence_colname].values,index=sequence_tbl[keycolname]).to_dict()
            else:
                seqlist = sequence_tbl[sequence_colname].tolist()
    else: #sequence_tbl is a list
        seqlist = list(sequence_tbl)
    return {"%s%d" % (keyname, (k+1)): v for k, v in enumerate(seqlist)}

def itoseq(seqint,kmer):
    nucleotides = {0:'A',1:'C',2:'G',3:'T'}
    seq = ""
    while(seqint > 0):
        seq = nucleotides[seqint & 3] + seq
        seqint >>= 2
    while len(seq) < kmer:
        seq = 'A' + seq
    return seq

'''
does not append 1, used for integer indexing
'''
def seqtoi(seq):
    nucleotides = {'A':0,'C':1,'G':2,'T':3}
    binrep = 0
    for i in range(0,len(seq)):
        binrep <<= 2
        binrep |= nucleotides[seq[i]]
    return binrep
