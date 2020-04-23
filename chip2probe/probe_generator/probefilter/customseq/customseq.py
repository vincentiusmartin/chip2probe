'''
This file contains functions used in making custom sequences.

Created on Aug 21, 2019

By: Vincentius Martin, Farica Zhuang
'''

import itertools
import pandas as pd
import random

def make_flank_products(input_sequence, start_site, end_site, mutable_flank_length = 1):
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

def make_weak_sites(df, escores, models, tf_source, ncustom):
    """
    Weaken sites in sequences.
    df: the dataframe containing sequences to weaken
    escores: dictionary of escore objects for the proteins
    models: dictionary of model objects for the proteins
    Return: a list of dictionary of seqeunces with sites weakened
    """
    weak_custom = []

    # here, assume sequence have 2 sites
    i = 0
    for index, row in df.iterrows():
        key = tf_source + "_" + row.key
        print("Working on %s" % key)
        seq = row.flank_left + row.wt + row.flank_right
        flank_len = len(row.flank_left)
        seq_len = len(seq)
        #print("  Making weak site for site1")
        w1 = make_weak_site(seq, escores[row.tf1], models[row.tf1],
                            which_site=0, mutable_flank_length=2,
                            left_outer=True) # left
        if seq == w1:
            print("could not find weak sites for the first site")
            continue
        #print("  Making weak site for site2")
        w2 = make_weak_site(seq, escores[row.tf2], models[row.tf2], 
                            which_site=1, mutable_flank_length=2, 
                            left_outer=False) # right
        if seq == w2:
            print("could not find weak sites for the second site")
            continue
        #print("  Making weak site for site1 and site2")
        w3 = make_weak_site(w1, escores[row.tf2], models[row.tf2], 
                            which_site=1, mutable_flank_length=2, 
                            left_outer=False) # both
        if w1 == w3 or seq == w3:
            print("could not find weak sites for both sites")
            continue
        
        weak_custom.append({"key": "%s_weak_s1" % key,
                            "flank_left": w1[0:flank_len],
                            "sequence": w1[flank_len:flank_len + seq_len],
                            "flank_right": w1[flank_len + seq_len:2 * flank_len + seq_len]})
        weak_custom.append({"key": "%s_weak_s2" % key,
                            "flank_left": w2[0:flank_len],
                            "sequence": w2[flank_len:flank_len + seq_len],
                            "flank_right": w2[flank_len + seq_len:2 * flank_len + seq_len]})
        weak_custom.append({"key": "%s_weak_s1_2" % key,
                            "flank_left": w3[0:flank_len],
                            "sequence": w3[flank_len:flank_len + seq_len],
                            "flank_right": w3[flank_len + seq_len:2 * flank_len + seq_len]})
        # if we have enough then break
        if len(weak_custom) >= ncustom * 3:
            break

    return weak_custom


def make_weak_site(input_sequence, pbmescore, imads, which_site, 
                   mutable_flank_length=2, left_outer=True,
                   elowcutoff=0.35, esignifcutoff=0.4, imads_upper_thres=0.5):
    """
    Create weak sites in the sequence.
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


def move_sites(df, direction="to_right", add_flank=False, 
               fixed_dist=False, cushion=5):
    """
    Move one binding site to the left or right, fixing the other binding site.

    direction: the direction to move a binding site to: to_right, to_left
    add_flank: bool whether to add flank after moving binding site
    """
    if direction != "to_right" and direction != "to_left":
        raise Exception("Direction can only be 'to_right' or 'to_left'.")

    dist_custom = []
    for row in df.itertuples():
        # get the sequence
        seq = row.wt
        threshold = row.ecutoff
        gap = row.egapthres

        # define flank
        flank = ""
        if add_flank:
            if direction == "to_right":
                flank = row.flank_left
            else:
                flank = row.flank_right

        # distance is number of base pairs between the end of site 1 and start of site 2
        dist = row.core2_start - row.core1_end - 1

        # get step size
        if direction == "to_right":
            step = abs(dist % 2 - 2)
        else:
            step = dist % 2 - 2

        # move binding site
        if fixed_dist:
            new_seq = move_site_fix_dist(seq, row.core1_end, row.core2_start,
                                         flank, step_to_right=step,
                                         add_flank=add_flank, cushion=cushion)
            fixed = "fixed"
        else:
            new_seq = move_single_bsite(seq, row.core1_end, row.core2_start,
                                        flank, step_to_right=step,
                                        add_flank=add_flank)
            fixed = ""

        # append the new custom sequence
        dist_custom.append({"key": "%s_%s_dist_%d" % (row.key, fixed, step),
                            "flank_left": row.flank_left, "sequence": new_seq,
                            "flank_right": row.flank_right})

        # update the step size
        if direction == "to_right":
            step += 2
        else:
            step -= 2

    return dist_custom


def move_site_fix_dist(input_sequence, start_site, end_site, flank_to_append, 
                       step_to_right=0, cushion=5, add_flank=True):
    """
    Move binding sites while fixing their distance.

    If step > 0, we're moving the binding sites to the right so add flank on the left
    Cushion is the minimum number of binding pairs at the ends
    """
    sequence = str(input_sequence)
    seqlen = len(input_sequence)

    # check that there is enough flank
    if add_flank and abs(step_to_right) > len(flank_to_append):
        raise Exception("step could not be larger than flank")
    # move binding site
    if step_to_right >= 0:
        if add_flank and end_site + step_to_right > seqlen:
            raise Exception("step is larger than allowed")
        movedseq = sequence[:seqlen - cushion - step_to_right] + sequence[seqlen - cushion:]
        if add_flank:
            newseq = flank_to_append[-step_to_right:] + movedseq
        else:
            newseq = movedseq
    else:
        step = abs(step_to_right)
        if add_flank and start_site + step_to_right < 0:
            raise Exception("step is less than allowed")
        movedseq = sequence[:cushion] + sequence[cushion + step_to_right:]
        if add_flank:
            newseq = movedseq + flank_to_append[:step_to_right]
        else:
            newseq = movedseq
    return newseq


def move_single_bsite(input_sequence, site1_end, site2_start, flank_to_append,
                      step_to_right=0, add_flank=True):
    """
    Shift the position of a binding site.

    Remove the middle sequence and append the flanking sequence of the same length removed.
    Take s1 as the bsite on the left and s2 as the bsite on the right
    If step_to_right is positive, move s1 to the right while fixing s2.
    If step_to_right is negative, move s2 to the left while fixing s1.
    Return: New sequence

    TODO: give an option to generate random sequence instead if flank is not known
    """
    # Get the original sequence and length of the sequence
    sequence = str(input_sequence)
    seqlen = len(input_sequence)

    # Check if there is enough flank provided
    if add_flank and abs(step_to_right) > len(flank_to_append):
        raise Exception("step could not be larger than flank")
    # check if the the shift is valid
    if step_to_right > 0 and site1_end + step_to_right > site2_start:
        raise Exception("step is larger than allowed")
    # check if the shift is valid
    if step_to_right < 0 and site2_start + step_to_right < site1_end:
        raise Exception("step is larger than allowed")

    # determine indices of the sequence to be removed
    mid = site1_end + (seqlen // 2)
    start = mid - (abs(step_to_right) // 2)
    end = mid + (abs(step_to_right) // 2)
    if abs(step_to_right) % 2 == 1:
        end += 1

    # get the sequence with a shifted binding site
    movedseq = sequence[:start] + sequence[end:]

    # add flank if needed
    if add_flank:
        if step_to_right >= 0:
            movedseq = flank_to_append[-(seqlen - len(movedseq)):] + movedseq
        else:
            movedseq = movedseq + flank_to_append[:seqlen - len(movedseq)]

    # return the new sequence
    return movedseq
