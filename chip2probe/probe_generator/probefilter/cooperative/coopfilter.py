"""
This file contains functions used in the filtering step of sequences.

Created on Aug 8, 2019

Authors: Vincentius Martin, Farica Zhuang
"""
from chip2probe.probe_generator.probefilter import sitespredict
from sitespredict.sitesplotter import SitesPlotter
from sitespredict.imads import iMADS
from sitespredict.sequence import Sequence
from sitespredict.kompas import Kompas
from cooperative import coopfilter
from chip2probe.util import bio as bio

from Bio.Seq import Seq
import itertools
import pandas as pd

def mutate_escore_seq_at_junc(seq, proteins, pos, comb_num, mutate_cutoff,
                              escores):
    """Return a mutated sequence that is nonspecific for all proteins."""
    if len(seq) != 8:
        raise Exception("sequence must be at length of 8")

    nucleotides = ["A", "C", "G", "T"]
    all_muts = {}

    # get every mutation for the sequence at the given position
    for comb in itertools.product(nucleotides, repeat=comb_num):
        seqlist = list(seq)
        # mutate the sequence
        seqlist[pos] = comb[0]
        if len(comb) > 1:
            seqlist[pos + 1] = comb[1]
        mutseq = "".join(seqlist)
        all_muts[mutseq] = {}
        # get esocre of this mutation for each protein
        for protein in proteins:
            # get the escore for the mutated sequence
            all_muts[mutseq][protein] = escores[protein].predict_sequence(mutseq).predictions[0]['score']

    # find a mutated sequence that is non-specific for all proteins
    min_sum = float("inf")
    min_seq = ""

    # iterate through mutated sequences
    for seq in all_muts:
        # get escores
        if all(i < mutate_cutoff for i in all_muts[seq].values()):
            if sum(all_muts[seq].values()) < min_sum:
                min_sum = sum(all_muts[seq].values())
                min_seq = seq

    return min_seq


def clean_junction_rev_seq(ori_seq, proteins, mutate_cutoff,
                           escores, max_mutate_count=2):
    """
    Mutate sequence at free end.

    Return: New sequence mutated at free end, reverse complement primer+new sequence
    """
    # get the first 8mer of the original sequence that we want to mutate
    seq = ori_seq
    curr_count = 1
    while seq == ori_seq and curr_count <= max_mutate_count:
        short_seq = seq[:8]
        # get the new sequence with mutated free end
        seq = mutate_escore_seq_at_junc(seq=short_seq,
                                        proteins=proteins,
                                        pos=curr_count - 1,
                                        comb_num=curr_count,
                                        mutate_cutoff=mutate_cutoff,
                                        escores=escores) + \
              seq[8:]
        curr_count += 1
    return seq


def clean_junction_seq(ori_seq, proteins, mutate_cutoff,
                       escores, max_mutate_count=2):
    """
    Mutate sequence at junction.

    count: max number of nucleotides to mutate
    Return: New sequence mutated at junction, new sequence+primer
    """
    seq = ori_seq
    curr_count = 1
    while seq == ori_seq and curr_count <= max_mutate_count:
        short_seq = seq[len(seq) - 8: len(seq)]
        # get the new sequence with the mutated 8mer
        seq = seq[:len(seq) - 8] + \
              mutate_escore_seq_at_junc(seq=short_seq,
                                        proteins=proteins,
                                        pos=len(short_seq) - curr_count,
                                        comb_num=curr_count,
                                        mutate_cutoff=mutate_cutoff,
                                        escores=escores)
        curr_count += 1
    return seq


def clean_junctions(seqlst, proteins, escores, models, mutate_cutoff=0.38,
                    mutate_gap=0, primer="GTCTTGATTCGCTTGACGCTGCTG",
                    max_mutate_count=2):
    """
    Check if junctions of the sequence are clean when combined with primer.

    Return: True if junctions or clean or if mutations are successful.
            False otherwise
    """
    # loop through each sequence in seqlst
    for i in range(len(seqlst)):
        seq = seqlst[i]
        # if there is not sequence, don't do anything
        if seq == "":
            continue
        es_preds = {}
        model_preds = {}
        es_preds_primer = {}
        model_preds_primer = {}
        es_preds_primer = {}
        model_preds_primer = {}
        d = {"sequence": [seq]}
        df = pd.DataFrame.from_dict(d)
        d_primer = {"sequence": [seq + primer]}
        df_primer = pd.DataFrame.from_dict(d_primer)
        # check for forward orientation
        for protein in proteins:
            # predict sequence without primer
            es_preds[protein] = escores[protein].predict_sequences(df)["sequence1"]
            model_preds[protein] = models[protein].predict_sequences(df)["sequence1"]
            # predict sequence with primer
            es_preds_primer[protein] = escores[protein].predict_sequences(df_primer)["sequence1"]
            model_preds_primer[protein] = models[protein].predict_sequences(df_primer)["sequence1"]
        seq_bs = Sequence(es_preds, model_preds, proteins, escores, mutate_cutoff, mutate_gap)
        seq_primer_bs = Sequence(es_preds_primer, model_preds_primer, proteins,
                                 escores, mutate_cutoff, mutate_gap)
        # theres a sequence with less bsites with primer. why?
        if seq_bs.site_count_all() > seq_primer_bs.site_count_all():
            print("Number of bsites with primer should not be less than without")
            seqlst[i] = ""
        if seq_bs.site_count_all() < seq_primer_bs.site_count_all():
            new_seq = clean_junction_seq(ori_seq=seq,
                                         proteins=proteins,
                                         mutate_cutoff=mutate_cutoff,
                                         escores=escores,
                                         max_mutate_count=2)
            if new_seq == seq:
                seqlst[i] = ""

    # check reverse complement
    for i in range(len(seqlst)):
        # get the sequence
        seq = seqlst[i]
       # get the reverse complement of the primer
        primer_obj = Seq(primer)
        rev_primer = str(primer_obj.reverse_complement())
        # if there is not sequence, don't do anything
        if seq == "":
            continue

        es_preds = {}
        model_preds = {}
        es_preds_primer = {}
        model_preds_primer = {}

        d = {"sequence": [seq]}
        df = pd.DataFrame.from_dict(d)
        d_primer = {"sequence": [rev_primer + seq]}
        df_primer = pd.DataFrame.from_dict(d_primer)

        # check for reverse complement
        for protein in proteins:
            # predict sequence without primer
            es_preds[protein] = escores[protein].predict_sequences(df)["sequence1"]
            model_preds[protein] = models[protein].predict_sequences(df)["sequence1"]
            # predict reverse complement of sequence with primer
            es_preds_primer[protein] = escores[protein].predict_sequences(df_primer)["sequence1"]
            model_preds_primer[protein] = models[protein].predict_sequences(df_primer)["sequence1"]
        seq_bs = Sequence(es_preds, model_preds, proteins, escores, mutate_cutoff, mutate_gap)
        seq_primer_bs = Sequence(es_preds_primer, model_preds_primer, proteins,
                                 escores, mutate_cutoff, mutate_gap)
        if seq_bs.site_count_all() > seq_primer_bs.site_count_all():
            print("Number of bsites with primer should not be less than without")
            seqlst[i] = ""
        if seq_bs.site_count_all() < seq_primer_bs.site_count_all():
            new_seq = clean_junction_rev_seq(ori_seq=seq,
                                             proteins=proteins,
                                             mutate_cutoff=mutate_cutoff,
                                             escores=escores,
                                             max_mutate_count=2)
            if new_seq == seq:
                seqlst[i] = ""

            else:
                seqlst[i] = seqlst[i] + primer

    # if wild type is empty, return false
    if seqlst[0] == "":
        return [], False

    # otherwise, return the new seqlist
    return seqlst, True


def check_sequence(seq, sites_list, pbmescores,
                   escore_cutoff=0.4, escore_gap=0,
                   get_complete_mutated=True):
    """
    Check if given sites are within specific 8mers in the sequence.

    seq: sequence to be checked
    sites_list: list of regions expected to be significant in sorted order
    pbmescore: pbmescore object
    escore_cutoff: escore_cutoff for significance
    escore_gap: number of gaps in significant escores allowed
    """
    # if the sequence is an empty string, return True by default
    if seq == "":
        if get_complete_mutated == False:
            return True
        return False
    # get the list of regions (represented as dictionaries) in the sequence
    # that are specific/significant in sorted order
    especific = []
    for escore in pbmescores:
        especific += escore.get_escores_specific(seq, escore_cutoff, escore_gap)
    # if number of specific regions expected and assumed are different, return false
    if len(sites_list) != len(especific):
        return False
    # for each significant region, check if if expected significant site falls inside
    for i in range(len(especific)):
        e_start = especific[i]["startpos"]
        e_end = especific[i]["startpos"] + especific[i]["escorelength"]
        core_start = sites_list[i][0]
        core_end = sites_list[i][1]
        # if not, return false
        if not (e_start <= core_end and core_start <= e_end):
            return False
    return True


def check_all_seqs(wt, m1, m2, m3, sites_dict, pbmescores, escore_cutoff=0.4,
                   escore_gap=0, get_complete_mutated=True):
    """
    Return True only if all conditions are met.

    Conditions:
    wt has two specific sites
    m1 has non-specific site at position 1
    m1 has specific site at position 2
    m2 has specific site at position 1
    m1 has non-specific site at position 2
    m3 has non specific sites at positions 1 and 2

    get_complete_mutated: if True then we require all 4 sequences
                          if False then we allow m1 m2 m3 to be empty
    """
    # check if all sequences are different
    if get_complete_mutated and len(set([wt, m1, m2, m3])) != 4:
        return False
    # initialize lists of the start and end positions of the two cores
    s1 = [sites_dict["core_start_1"], sites_dict["core_end_1"]]
    s2 = [sites_dict["core_start_2"], sites_dict["core_end_2"]]
    # initialize a list of the pbmescores for the two cores in order
    if len(pbmescores) > 1:
        if s1[0] < s2[0]:
            pbmescores = [pbmescores[sites_dict["protein_1"]], pbmescores[sites_dict["protein_2"]]]
        else:
            pbmescores = [pbmescores[sites_dict["protein_2"]], pbmescores[sites_dict["protein_1"]]]
    else:
        pbmescores = list(pbmescores.values())
    # for wild type, both sites are specific
    wt_cond = check_sequence(wt, [s1, s2], pbmescores, escore_cutoff,
                             escore_gap, get_complete_mutated=False)
    # for m1, only s2 is specific
    m1_cond = check_sequence(m1, [s2], pbmescores, escore_cutoff,
                             escore_gap, get_complete_mutated=get_complete_mutated)
    # for m2, only s1 is specific
    m2_cond = check_sequence(m2, [s1], pbmescores, escore_cutoff,
                             escore_gap, get_complete_mutated=get_complete_mutated)
    # for m3, both sites are non-specific
    m3_cond = check_sequence(m3, [], pbmescores, escore_cutoff,
                             escore_gap, get_complete_mutated=get_complete_mutated)

    return (wt_cond and m1_cond and m2_cond and m3_cond)


def get_filtered_probes(seqdf, escores, models, mutate_cutoff, mutate_gap,
                        egaps, thresholds, proteins, colors,
                        generate_plots=False, spcomb=[(0, 0)], analysis_path="",
                        mode="custom", predict_flanks=True, flank_len=0,
                        key_colname="key",
                        show_model_flanks=False, get_complete_mutated=True,
                        primer="", max_mutate_count=2):
    """Get the filtered probes with m1,m2,m3 for each sequence in the given df."""
    filtered_probes = []
    # iterate through each site num and peak len combination
    for comb in spcomb:
        # get escore and model predictions for each protein
        es_preds = {}
        esplots = {}
        model_preds = {}
        model_plots = {}
        sitenum = comb[0]
        peaklen = comb[1]

        # get rows with the current sitenum and peaklen if specified
        if sitenum != 0 and peaklen != 0:
            df = seqdf.loc[(seqdf["sites_in_peak"] == sitenum) & (seqdf["peaklen"] == peaklen)]
        # otherwise use all rows
        else:
            df = seqdf
        # initialize escore and model objects for each protein
        for protein in proteins:
            protein_num = proteins.index(protein)
            es_preds[protein] = escores[protein].predict_sequences(df, key_colname=key_colname)
            esplots[protein] = escores[protein].make_plot_data(es_preds[protein], color=colors[protein_num][0])

            model_preds[protein] = models[protein].predict_sequences(df,
                                                                     key_colname=key_colname,
                                                                     predict_flanks=predict_flanks,
                                                                     flank_len=flank_len)
            model_plots[protein] = models[protein].make_plot_data(model_preds[protein],
                                                                  color=colors[protein_num][1],
                                                                  show_model_flanks=show_model_flanks)

        # Generate plots
        if generate_plots:
            sp = SitesPlotter()
            # if need to plot, uncomment this
            sp.plot_seq_combine([esplots, model_plots],
                                filepath="%s/sitesplot_d%d_p%d.pdf" %
                                (analysis_path, sitenum, peaklen))

        # get filtered sequences
        filtered_seqs = {}
        flanks = {}
        print("Site filtering...")
        print("Number of sites before mutating:", len(es_preds[proteins[0]]))

        # get sequences with 2 significant binding sites
        sites_mutated = 0
        sites_removed = 0
        failed_mutations = 0
        for key in es_preds[proteins[0]]:
            curr_es_preds = {}
            curr_model_preds = {}
            for protein in proteins:
                curr_es_preds[protein] = es_preds[protein][key]
                curr_model_preds[protein] = model_preds[protein][key]
            #print(key,"asd",curr_model_preds["ets1"])
            bs = Sequence(curr_es_preds, curr_model_preds, proteins=proteins,
                          escore_cutoff=mutate_cutoff, escore_gap=mutate_gap,
                          pbmescores=escores)
            ### print(key, bs.is_valid())
            if bs.is_valid():
                filtered_seqs[key] = bs
        # TODO: move all print statements to a log file
        # print("Number of sites mutated:", sites_mutated)
        # print("Number of failed mutations:", failed_mutations)
        # print("Number of sites removed:", sites_removed)
        print("Number of sites after filtering:", len(filtered_seqs))

        print("Creating m1,m2,m3 sequences...")
        # for each of the filtered sequence, create m1,m2,m3 sequences
        seqdict = {}
        funcdict = {}
        for key in filtered_seqs:
            # Visualization part
            seqdict["%s-wt" % key] = filtered_seqs[key].sequence
            # current binding site object
            bs = filtered_seqs[key]
            # get m1,m2,m3 for each wt
            for idx, mut in enumerate([[0], [1], [0, 1]]):
                # here we mutate on the first, second, and both sites
                # mut is the index of the site to abolish
                to_remove = bs.remove_pos(mut)
                mutseq = bs.abolish_sites(to_remove, mode="to_eliminate",
                                          escore_threshold=mutate_cutoff)
                seqdict["%s-m%d" % (key, idx + 1)] = mutseq.sequence
                funcdict["%s-m%d" % (key, idx + 1)] = mutseq.plot_functions

            # get sequences that pass given escore gap and threshold combination
            for e in list(itertools.product(egaps, thresholds)):
                egapthres = e[0]
                ecutoff = e[1]

                # check that wt, m1, m2, m3 are valid
                if coopfilter.check_all_seqs(seqdict["%s-wt" % key],
                                             seqdict["%s-m1" % key],
                                             seqdict["%s-m2" % key],
                                             seqdict["%s-m3" % key],
                                             filtered_seqs[key].get_sites_dict(),
                                             escores,
                                             escore_cutoff=ecutoff,
                                             escore_gap=egapthres,
                                             get_complete_mutated=get_complete_mutated):
                    bsites_dict = filtered_seqs[key].get_sites_dict()
                    lst = [seqdict["%s-wt" % key], seqdict["%s-m1" % key], seqdict["%s-m2" % key],
                           seqdict["%s-m3" % key]]
                    lst, successful = clean_junctions(seqlst=lst,
                                                      proteins=proteins,
                                                      escores=escores,
                                                      models=models,
                                                      mutate_cutoff=mutate_cutoff,
                                                      mutate_gap=mutate_gap,
                                                      primer="GTCTTGATTCGCTTGACGCTGCTG",
                                                      max_mutate_count=max_mutate_count)
                    if successful:
                        # replace seqdict with the new sequences
                        seqdict["%s-wt" % key] = lst[0]
                        seqdict["%s-m1" % key] = lst[1]
                        seqdict["%s-m2" % key] = lst[2]
                        seqdict["%s-m3" % key] = lst[3]
                        filtered_probes.append({"key": key,
                                                "wt": seqdict["%s-wt" % key],
                                                "m1": seqdict["%s-m1" % key],
                                                "m2": seqdict["%s-m2" % key],
                                                "m3": seqdict["%s-m3" % key],
                                                "tf1": bsites_dict["protein_1"],
                                                "tf2": bsites_dict["protein_2"],
                                                "core1_start": bsites_dict["core_start_1"],
                                                "core1_mid": bsites_dict["core_mid_1"],
                                                "core1_end": bsites_dict["core_end_1"],
                                                "core1_pref": bsites_dict["score_1"],
                                                "core2_start": bsites_dict["core_start_2"],
                                                "core2_mid": bsites_dict["core_mid_2"],
                                                "core2_end": bsites_dict["core_end_2"],
                                                "core2_pref": bsites_dict["score_2"],
                                                "ecutoff": ecutoff,
                                                "egapthres": egapthres,
                                                "distance": filtered_seqs[key].get_sites_dist(),
                                                "sites_in_peak": sitenum,
                                                "peak_length": peaklen
                                                })
                        break # the sequence passes the filtering check, so stop

        # generate plots of wt, m1, m2, m3
        if generate_plots:
            filtered_es_preds = {}
            filtered_esplots = {}
            filtered_model_preds = {}
            filtered_model_plots = {}
            for protein in proteins:
                protein_num = proteins.index(protein)
                filtered_es_preds[protein] = escores[protein].predict_sequences(seqdict, key_colname="key")
                filtered_esplots[protein] = escores[protein].make_plot_data(filtered_es_preds[protein], color=colors[protein_num][0])

                filtered_model_preds[protein] = models[protein].predict_sequences(seqdict,
                                                                                  key_colname="key",
                                                                                  predict_flanks=predict_flanks)
                filtered_model_plots[protein] = models[protein].make_plot_data(filtered_model_preds[protein],
                                                                               color=colors[protein_num][1],
                                                                               show_model_flanks=show_model_flanks)
                sp.plot_seq_combine([filtered_esplots, filtered_model_plots],
                                   filepath="%splot_%s_d%d_p%d.pdf" % (analysis_path, mode, sitenum, peaklen))

    return filtered_probes
