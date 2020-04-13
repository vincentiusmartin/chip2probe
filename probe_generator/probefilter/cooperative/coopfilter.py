"""
Created on Aug 8, 2019.

Authors: Vincentius Martin, Farica Zhuang
"""


def are_specific_8mers_within_coords(seq, coord_list, pbmescore,
                                     escore_cutoff=0.4, escore_gap=0):
    """
    Determine if given regions have specific 8mers.

    seq: sequence to be checked
    coord_list: list of regions expected to be significant in sorted order
    pbmescore: pbmescore object
    escore_cutoff: escore_cutoff for significance
    escore_gap: number of gaps in significant escores allowed
    """
    # if the sequence is an empty string, return True by default
    if seq == "":
        return True
    # get the list of regions (represented as dictionaries) in the sequence
    # that are specific/significant in sorted order
    especific = pbmescore.get_escores_specific(seq, escore_cutoff, escore_gap)
    # if number of specific regions expected and assumed are different, return false
    if len(coord_list) != len(especific):
        return False
    # for each significant region, check if if expected significant site falls inside
    for i in range(len(especific)):
        e_start = especific[i]["startpos"]
        e_end = especific[i]["startpos"] + especific[i]["escorelength"]
        core_start = coord_list[i][0]
        core_end = coord_list[i][1]
        # if not, return false
        if not (e_start <= core_end and core_start <= e_end):
            return False
    return True

def filter_coopseq(wt, m1, m2, m3, sites_dict, pbmescore, escore_cutoff=0.4, escore_gap = 0):
    """
    return True only if all conditions are met:
    wt has two specific sites
    m1 has non-specific site at position 1
    m1 has specific site at position 2
    m2 has specific site at position 1
    m1 has non-specific site at position 2
    m3 has non specific sites at positions 1 and 2
    """

    # check if all sequences are different
    if len(set([wt, m1, m2, m3])) != 4:
        return False
    # initialize lists of the start and end positions of the two cores
    s1 = [sites_dict["core_start_1"], sites_dict["core_end_1"]]
    s2 = [sites_dict["core_start_2"], sites_dict["core_end_2"]]
    # for wild type, both sites are specific
    wt_cond = are_specific_8mers_within_coords(wt, [s1, s2], pbmescore,
                                               escore_cutoff, escore_gap)
    # for m1, only s2 is specific
    m1_cond = are_specific_8mers_within_coords(m1, [s2], pbmescore,
                                               escore_cutoff, escore_gap)
    # for m2, only s1 is specific
    m2_cond = are_specific_8mers_within_coords(m2, [s1], pbmescore,
                                               escore_cutoff, escore_gap)
    # for m3, both sites are non-specific
    m3_cond = are_specific_8mers_within_coords(m3, [], pbmescore,
                                               escore_cutoff, escore_gap)

    return (wt_cond and m1_cond and m2_cond and m3_cond)
