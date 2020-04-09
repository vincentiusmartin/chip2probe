'''
Created on Aug 21, 2019

By: Vincentius Martin, Farica Zhuang
'''

import itertools

import random

def make_flank_products(input_sequence, start_site, end_site, mutable_flank_length = 1):
    '''Get '''
    #seuence = str(input_sequence)
    flank_pos = [*range(start_site - mutable_flank_length, start_site), *range(end_site, end_site + mutable_flank_length)]

    bases = ['A','C','G','T']
    flankcomb = list(itertools.product(bases, repeat = len(flank_pos)))
    all = {}
    for i in range(len(flankcomb)):
        curseq = list(input_sequence)
        for pos,base in zip(flank_pos,flankcomb[i]):
            curseq[pos] = base
        #print(flankcomb[i],curseq)
        all["".join(flankcomb[i])] = "".join(curseq)

    return all


def make_weak_site(input_sequence, pbmescore, imads, which_site, mutable_flank_length = 2, left_outer = True,
                   elowcutoff=0.35, esignifcutoff=0.4, imads_upper_thres = 0.5):
    """
    which_site = site index
    looseness_count = how many escore can b below cutofff
    TODO enable left_outer to be fully random as well
    """
    prevseqs = []

    ims_pred = imads.predict_sequence(input_sequence)
    if which_site >= len(ims_pred):
        print("There is no site %d" % which_site)
        return str(input_sequence)
    ims_init = ims_pred[which_site]
    if ims_init["score"] < imads_upper_thres:
        print("The site is already a weak site since it is below imads upper threshold")
        return str(input_sequence)
    site_start = ims_init['core_start']
    site_end = site_start + ims_init['core_width']

    # we start by mutating the outer first, if left_outer is true then start by mutating the left position
    flank_left = list(range(site_start - mutable_flank_length, site_start))
    flank_right = list(range(site_end, site_end + mutable_flank_length))
    flank_pos = [*flank_left, *flank_right]

    # make ordering
    wt_flank = [input_sequence[p] for p in flank_pos]
    bases = ['A','C','G','T']
    mid = len(wt_flank)  //  2
    flankcomb = []
    for i in range(len(wt_flank)):
        flen = i + 1
        # still not effective as we call this every time, fix later
        combs = list(itertools.product(bases, repeat = flen))
        random.shuffle(combs)
        cur = []
        if left_outer:
            toleft = 0 if mid - flen < 0 else mid - flen
            toright = mid if i - mid < 0 else flen
        else:
            toleft = mid if mid - flen >= 0 else len(wt_flank) - flen
            toright = mid + flen if i - mid < 0 else len(wt_flank)
        for comb in combs:
            cur_flank = list(wt_flank)
            cur_flank[toleft:toright] = comb
            cf = "".join(cur_flank)
            if cf not in flankcomb:
                flankcomb.append(cf)

    flag = False
    i = 0
    curseq = ""
    while not flag and i < len(flankcomb):
        # make the sequence
        curseq = list(input_sequence)
        for pos,base in zip(flank_pos, flankcomb[i]) :
            curseq[pos] = base
            #print(flankcomb[i],curseq)
        curseq = "".join(curseq)
        ims_preds = imads.predict_sequence(curseq)
        if which_site < len(ims_preds):
            ims = imads.predict_sequence(curseq)[which_site]

            # check if we get the same sequence
            if ims and ims['core_start'] == site_start and ims['core_start'] + ims['core_width'] == site_end:
                epreds = pbmescore.predict_sequence(curseq)
                epreds_in_core = list(filter(lambda epred:
                                                epred['position'] >= site_start and
                                                epred['position'] < site_end,
                                            epreds.predictions))
                e_within_range = [True if e['score'] >= elowcutoff and e['score'] < esignifcutoff else False for e in epreds_in_core]
                # get 2 consecutive e-scores  within range
                pass_erange = any(p1 and p2 for p1,p2 in zip(e_within_range, e_within_range[1:]))
                flag = ims['score'] < imads_upper_thres and ims['score'] > imads.imads_threshold  and pass_erange
        i += 1
    if flag:
        return curseq
    else:
        print("Could not find any weak sites with the given criteria")
        return str(input_sequence)

def move_single_bsite(input_sequence, site1_end, site2_start, flank_to_append, step_to_right = 0, add_flank = True):
    '''
    Shift the position of a binding site.
    Remove the middle sequence and append the flanking sequence of the same length removed.
    Take s1 as the bsite on the left and s2 as the bsite on the right
    If step_to_right is positive, move s1 to the right while fixing s2.
    If step_to_right is negative, move s2 to the left while fixing s1.
    Return: New sequence

    TODO: give an option to generate random sequence instead if flank is not known
    '''

    # Get the original sequence and length of the sequence
    sequence = str(input_sequence)
    seqlen = len(input_sequence)

    # Check if there is enough flank provided
    if add_flank and abs(step_to_right) > len(flank_to_append):
        raise Exception("step could not be larger than flank")
     # check if the the shift is valid
    if step_to_right > 0 and site1_end + step > site2_start:
        raise Exception("step is larger than allowed")  
    # check if the shift is valid
    if step_to_right < 0 and site2_start + step < site1_end:
        raise Exception("step is larger than allowed")

    mid = site1_end+(seqlen//2)
    start = mid-(step//2)
    end = mid+(step//2) 
    if step%2==1:
        end += 1
    # get the sequence
 with a shifted binding site
    movedseq = sequence[:start] + sequence[end:]

    # add flank if needed
    if add_flank:
        if step_to_right >= 0:
            movedseq = flank_to_append[-(seqlen-len(movedseq)):] + movedseq
        else:
            movedseq = movedseq + flank_to_append[:seqlen-len(movedseq)]

    # return the new sequence
    return movedseq
