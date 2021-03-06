'''
Created on Aug 8, 2019

@author: Farica Zhuang, Vincentius Martin
'''

def are_specific_8mers_within_coords(seq, coord_list, pbmescores, 
                                     escore_cutoff=0.4, escore_gap = 0, 
                                     left=""):
    """
    Check if the 8mers within the cores are specific."""
    especific = []
    for escore in pbmescore:
        especific += pbmescore.get_escores_specific(seq, escore_cutoff, escore_gap)
    if len(coord_list) != len(especific):
        return False
    # for now assume both are coord_list and especific results are sorted
    for i in range(len(especific)):
        e_start = especific[i]["startpos"]
        e_end = especific[i]["startpos"] + especific[i]["escorelength"]
        core_start = coord_list[i][0]
        core_end = coord_list[i][1]
        if not (e_start <= core_end and core_start <= e_end):
            return False
    return True

def filter_coopseq(wt, m1, m2, m3, sites_dict, pbmescores, escore_cutoff=0.4, escore_gap = 0):
    """
    return True only if some conditions are met
    this assume that: 
    m1 is mutation at position 1
    m2 is mutation at position 2
    m3 is mutation at position 1 and 2
    """
    # check if all sequences are different
    if len(set([wt,m1,m2,m3])) != 4:
        return False
    s1 = [sites_dict["core_start_1"],sites_dict["core_end_1"]]
    s2 = [sites_dict["core_start_2"],sites_dict["core_end_2"]]
    if s1[0] < s2[0]:
            left = s1["protein"]
            right = s2["protein"]
        else:
            left = s2["protein"]
            right = s2["protein"]
    if len(pbmescores) > 1:
        pbmescores = [pbmescores[left], pbmescores[right]]
    else:
        pbmesocres = [pbmescores[left], pbmescores[left]]
    wt_cond = are_specific_8mers_within_coords(wt, [s1,s2], pbmescores, escore_cutoff, escore_gap)
    m1_cond = are_specific_8mers_within_coords(m1, [s2], pbmescores, escore_cutoff, escore_gap)
    m2_cond = are_specific_8mers_within_coords(m2, [s1], pbmescores, escore_cutoff, escore_gap)
    m3_cond = are_specific_8mers_within_coords(m3, [], pbmescores, escore_cutoff, escore_gap)
    return (wt_cond and m1_cond  and m2_cond and m3_cond)
    