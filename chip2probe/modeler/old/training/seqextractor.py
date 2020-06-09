'''
Created on Nov 1, 2019

@author: vincentiusmartin
'''

import itertools

# ========== Sequence extractor =============
def extract_kmer_ratio(seq,kmer):
    nucleotides = ['A','C','G','T']

    kmer_count = {}
    ratio = {}
    perm = ["".join(p) for p in itertools.product(nucleotides, repeat=kmer)]
    total = 0
    for p in perm:
        kmer_count[p] = 0
    for i in range(0,len(seq)+1-kmer):
        kmer_count[seq[i:i+kmer]] += 1
        total += 1
    for p in sorted(perm):
        if total > 0:
            ratio[p] = float(kmer_count[p])/total
        else:
            # can be zero if there is no linker
            ratio[p] = 0
    return ratio

def extract_positional(seq, maxk = 2, label="seq", orientation="right"):
    '''
    orientation: if right, then start from 0 to the right, else start from
    len(seq)-1 to the left
    '''
    if orientation == "right":
        iterseq = str(seq)
    else:
        iterseq = str(seq[::-1])
    nucleotides = ['A','C','G','T']
    features = {}
    for k in range(1,maxk+1):
        perm = ["".join(p) for p in itertools.product(nucleotides, repeat=k)]
        for i in range(0,len(iterseq)+1-k):
            m = 1 if orientation == "right" else -1
            for kmer in perm:
                if orientation == "right":
                    seqcmp = iterseq[i:i+k]
                    idx = m * i + k - 1
                else:
                    seqcmp = iterseq[i:i+k][::-1]
                    idx = m * i - k + 1

                if seqcmp == kmer:
                    features["pos_%s_%d_%s" % (label,idx,kmer)] = 1
                else:
                    features["pos_%s_%d_%s" % (label,idx,kmer)] = 0
    return features

def extract_positional_features(seq, bpos1, bpos2, span_out, span_in):
    '''
    If the motif length is even, then it is important to insert bpos2 + 1 so
    we get the same span_in, span_out length on bpos1 and bpos2
    '''
    nucleotides = ['A','C','G','T']

    pos1 = bpos1 - 1 #
    pos2 = bpos2 - 1

    # start position
    b1_left = seq[pos1-span_out:pos1+1]
    b1_right = seq[pos1:pos1+span_in]

    b2_left = seq[pos2-span_in:pos2+1]
    b2_right = seq[pos2:pos2+span_out+1]

    feature1_l = extract_positional(b1_left,2,label="left",orientation="left")
    feature1_r = extract_positional(b1_right,2,label="left",orientation="right")
    feature2_l = extract_positional(b2_left,2,label="right",orientation="left")
    feature2_r = extract_positional(b2_right,2,label="right",orientation="right")

    return {**feature1_l, **feature1_r, **feature2_l, **feature2_r}
