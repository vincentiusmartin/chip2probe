import pandas as pd

# all permutations are already reverse-deleted
# all sequences are represented in binary


def revcompstr(seq):
    rev = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([rev[base] for base in reversed(seq)])


def get_seqdict(sequence_tbl, sequence_colname = "sequence", keyname="sequence"):
    """
    Generates dictionary representation from sequence table or list
    :param sequence_tbl: 
    :param sequence_colname: 
    :return: 
    """
    if not isinstance(sequence_tbl, pd.DataFrame) and not isinstance(sequence_tbl, list):
        raise Exception('sequence_tbl must be either pandas data frame or a list')
    if type(sequence_tbl) == pd.DataFrame:
        if not sequence_colname in sequence_tbl:
            raise Exception('Could not find column %s in the data frame' % sequence_colname)
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
