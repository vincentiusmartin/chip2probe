'''
Created on July 3, 2020

@author: vincentiusmartin, Farica Zhuang
@description: all general function for the cooperative projects
'''

import math


def move_single_site(sequence, site1_end, site2_start,
                    step, flank_append = "", patch=False):
    """
    Make sites closer by shifting the position of a binding site.

    Remove the middle sequence and append the flanking sequence of the same length removed.
    Take s1 as the bsite on the left and s2 as the bsite on the right
    If step is positive, move s1 to the right while fixing s2.
    If step is negative, move s2 to the left while fixing s1.

    Args:
        sequence: input sequence
        site1_end: the end position of site 1
        site2_start: the start position of site 2
        step: move the site to the left if positive, else to the right
        flank_append: flank sequence to append; if step is positive append left, else right.
        patch: if True, use the sequence removed to append

    Return:
        New sequence with site1 or site2 moved
    Example:
        seq = "CAGAGTAGGAAGCAGCTGTTCTAACTTCCGTCTTCG"
        move_single_site(sequence, 15, 21, 4,patch=True)

    TODO:
        give an option to generate random flank_append
    """
    # Get the original sequence and length of the sequence
    seqlen = len(sequence)

    # check if the the shift is valid
    if (step > 0 and site1_end + step > site2_start) or (step < 0 and site2_start + step < site1_end):
        raise Exception("step is larger than allowed")
    if len(flank_append) > 0 and patch:
        raise Exception("Should only specify either flank or patch")

    # determine indices of the sequence to be removed
    mid = (site1_end + site2_start) // 2
    mid_start = mid - (math.floor(abs(step) / 2))
    mid_end = mid_start + abs(step)
    if patch:
        flank = sequence[mid_start:mid_end]
    else:
        flank = flank_append

    if step > 0:
        movedseq = flank + sequence[:mid_start] + sequence[mid_end:]
    else:
        movedseq = sequence[:mid_start] + sequence[mid_end:] + flank

    return movedseq
