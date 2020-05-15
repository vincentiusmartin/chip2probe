"""
This file contains ProbeFilter class used to generate clean probes.

Authors: Farica Zhuang and Vincentius Martin
Created: April 21, 2020
"""
# these paths need to be handled better
import sys
import pandas as pd
import os
#sys.path.append(os.path.realpath("sitespredict/libsvm/python"))
from chip2probe.probe_generator.probefilter import sitespredict
from sitespredict.pbmescore import PBMEscore
from sitespredict.imadsmodel import iMADSModel
from sitespredict.imads import iMADS
from sitespredict.kompas import Kompas
from cooperative import coopfilter
from customseq import customseq
from functools import reduce


class ProbeFilter:
    """This class generate clean noncustom and custom probes."""

    def __init__(self, inseq_path, outpath="", **kwargs):
        """Initialize class variables."""
        for k, v in kwargs.items():
            setattr(self, k, v)
        if len(self.colors) == 0 or not self.colors[0]:
            self.colors = [('crimson', 'plum'), ('steelblue', 'lightblue')]
        self.inpath = inseq_path # "sites_all.tsv"
        self.analysis_path = outpath if outpath else os.getcwd()
        self.req_cols = ["key", "wt", "m1", "m2", "m3", "flank_left",
                         "flank_right", "core1_start", "core1_end",
                         "core1_mid", "core1_pref", "core2_start",
                         "core2_end", "core2_mid", "core2_pref",
                         "ecutoff", "egapthres", "distance",
                         "sites_in_peak", "peak_length",
                         "coordinate", "tf1", "tf2"]

    def _run_all(self):
        # initialize escore and model objects
        self.initialize_objects()

        # get filtered probes
        self.filter_probes()

        # # get custom probes
        # if len(self.customs) > 0:
        #     self.customize_probes()

        # # get negative controls
        # if self.num_neg_ctrl > 0:
        #     self.add_neg_ctrls()

        print("Done!!!!!")

    def initialize_objects(self):
        """Initialize escore and model objects."""
        escores = {}
        models = {}
        for tf in self.tf:
            escores[tf] = PBMEscore(self.escore_long_paths[tf])
            # initialize binding site prediction model used based on method
            if self.sitecall_mode == 'imads':
                # read imads prediction
                imads_models = [iMADSModel(modelpath, modelcore, self.imads[tf]["width"], [1, 2, 3])
                                for modelpath, modelcore in
                                zip(self.imads[tf]["model_paths"], self.imads[tf]["cores"])]
                models[tf] = iMADS(imads_models, self.imads[tf]["cutoff"])
            elif self.sitecall_mode == 'kompas':
                models[tf] = Kompas(protein=tf, threshold=self.mutate_cutoff,
                                         kmer_align_path=self.kompas[tf]["align_path"])
            # raise exception if model is not supported
            else:
                raise Exception("Supported methods are imads and kompas. Got {}".format(method))
        self.escores = escores
        self.models = models

    def filter_probes(self):
        """Filter and clean up probes so each sequence has only 2 significant sites."""
        ext = os.path.splitext(self.inpath)[1]
        sep = "\t" if ext == "tsv" else ","
        seqdf = pd.read_csv(self.inpath, sep=sep)

        # get all combinations of number of sites in peak and peak length as a list of tuples
        if "sites_in_peak" in seqdf:
            spcomb = sorted(list(set(zip(seqdf["sites_in_peak"], seqdf["peaklen"]))), key=lambda x: (x[0],x[1]))
        else:
            spcomb = [(0,0)] # default

        keycolname = "key" if "key" in seqdf else ""

        # python3 main.py /Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/chip2probe/modeling/mutlist.csv
        # for each combination, get filtered probes
        filtered_probes = coopfilter.get_filtered_probes(seqdf=seqdf,
                                                         escores=self.escores,
                                                         models=self.models,
                                                         mutate_cutoff=self.mutate_cutoff,
                                                         mutate_gap=self.mutate_gap,
                                                         egaps=self.escore_gaps,
                                                         thresholds=self.escore_thresholds,
                                                         proteins=self.tf,
                                                         colors=self.colors,
                                                         generate_plots=self.generate_plots,
                                                         spcomb=spcomb,
                                                         analysis_path=self.analysis_path,
                                                         mode="noncustom",
                                                         predict_flanks=True,
                                                         flank_len=self.flank_len,
                                                         key_colname=keycolname,
                                                         show_model_flanks=self.show_model_flanks,
                                                         get_complete_mutated=self.get_complete_mutated,
                                                         primer=self.primer,
                                                         max_mutate_count=self.max_mutate_count)

        # initialize a dataframe of filtered probes with wt, m1, m2, m3
        fp_df = pd.DataFrame(filtered_probes)

        df_init = seqdf[["key", "flank_left", "flank_right", "chr", "seqstart", "seqend"]] # will produce warning but it is fine
        df_init["coordinate"] = df_init.apply(lambda row: "%s:%d-%d" % (row["chr"],row["seqstart"],row["seqend"]), axis=1)
        df_init = df_init.drop(['chr', 'seqstart', 'seqend'], axis=1)
        df_out = pd.concat([fp_df.set_index('key'), df_init.set_index('key')], axis=1, join='inner').reset_index()

        # check if we have the required columns
        if df_out.columns.isin(self.req_cols).all():
            print("Saving result in %s/mutated_probes.tsv" % self.analysis_path)
            df_out.to_csv("%s/mutated_probes.tsv" % (self.analysis_path), sep="\t", index=False)
        else:
            print("No mutation rows found in %s" % sitepath)

    def customize_probes(self):
        """Customize clean probes."""
        infile = "%smutated_probes.tsv" % self.analysis_path
        df = pd.read_csv(infile, sep="\t")
        ret = []
        if 1 in self.customs:
            weak_sites = customseq.make_weak_sites(df, self.escores, self.models,
                                                   self.tf_source, self.ncustom)
            ret = ret + weak_sites
        if 2 in self.customs:
            moved_sites = customseq.move_sites(df, direction="to_right", add_flank=True,
                                               fixed_dist=False)
            ret = ret + moved_sites
        if 3 in self.customs:
            moved_sites = customseq.move_sites(df, direction="to_left", add_flank=True,
                                               fixed_dist=False)
            ret = ret + moved_sites
        if 4 in self.customs:
            moved_sites = customseq.move_sites(df, direction="to_right", add_flank=False,
                                               fixed_dist=False)
            ret = ret + moved_sites
        if 5 in self.customs:
            moved_sites = customseq.move_sites(df, direction="to_right", add_flank=True,
                                               fixed_dist=True)
            ret = ret + moved_sites
        if 6 in self.customs:
            moved_sites = customseq.move_sites(df, direction="to_left", add_flank=True,
                                               fixed_dist=True)
            ret = ret + moved_sites
        if 7 in self.customs:
            moved_sites = customseq.move_sites(df, direction="to_right", add_flank=False,
                                               fixed_dist=True)
            ret = ret + moved_sites


        #pickle.dump(weak_custom + dist_custom, open("test.pickle","wb"))
        #df_comb = pd.DataFrame(pickle.load(open("test.pickle","rb")))

        # turn the list of custom sequences to dataframe
        df_comb = pd.DataFrame(ret)

        # filter the sequences and obtain m1,m2,m3
        filtered_probes = coopfilter.get_filtered_probes(seqdf=df_comb,
                                                         escores=self.escores,
                                                         models=self.models,
                                                         mutate_cutoff=self.mutate_cutoff,
                                                         mutate_gap=self.mutate_gap,
                                                         egaps=self.egaps,
                                                         thresholds=self.thresholds,
                                                         proteins=self.proteins,
                                                         colors=self.colors,
                                                         generate_plots=self.generate_plots,
                                                         analysis_path=self.analysis_path,
                                                         mode="custom",
                                                         predict_flanks=True,
                                                         key_colname="key",
                                                         show_model_flanks=self.show_model_flanks,
                                                         get_complete_mutated=self.get_complete_mutated,
                                                         primer=self.primer,
                                                         max_mutate_count=self.max_mutate_count)

        # probably should check here if filtered_probes is empty
        fp_df = pd.DataFrame(filtered_probes)
        if fp_df.columns.isin(self.req_cols).all():
            print("Saving result in custom_probes_dist_and_weak.tsv")
            fp_df.to_csv("custom_probes_dist_and_weak.tsv", sep="\t", index=False)
        else:
            print("Failed to save custom probe file")

    def add_neg_ctrls(self):
        """Get negative controls"""
        ret = []
        # loop through the files provided
        for path in self.neg_ctrl_files:
            # get escore for runx for runx negative controls
            df = pd.read_csv(path)
            lst = self.get_neg_ctrls(df, self.num_neg_ctrl-len(ret))
            ret = ret + lst
            if len(ret) == self.num_neg_ctrl:
                break

        if len(ret) < self.num_neg_ctrl:
            for cell_line in self.dnase_paths:
                self.process_dnase(cell_line)

        ret_df = pd.DataFrame.from_records(list(ret))
        ret_df.to_csv("negative_controls.csv", index=None)

        print("Number of neg ctrls required:", self.num_neg_ctrl)
        print("Number of neg ctrls found:", len(ret_df))

    def get_neg_ctrls(self, df, num_seqs):
        """Get negative control seqs from given df."""
        es_preds = {}
        # get the escore prediction for each protein
        for protein in self.proteins:
            es_preds[protein] = self.escores[protein].predict_sequences(df,
                                                                        predict_flanks=False,
                                                                        sequence_colname="Sequence")
        # initialize lst to strong neg control sequences
        negs = set()
        # check each sequence in the df
        for key in es_preds[self.proteins[0]]:
            # check each escore in the sequence
            passed = True
            for protein in self.proteins:
                for pos in es_preds[protein][key].predictions:
                    # if any escore is at least cutoff or above,
                    # stop looping and indicated False
                    if pos["score"] >= self.neg_ctrl_thres:
                        passed = False
                        break
            # if the sequence passed the check, add to the list
            if passed:
                negs.add(es_preds[self.proteins[0]][key].sequence)
            # create a dataframe for runx negative controls
            # with no significant ets escore
            if len(negs) == num_seqs:
                break

        return list(negs)

    def process_dnase(self, cell_line):
        """Get negative control sequences from dnase seq."""
        # get dnase peaks
        dnase_peaks = self.get_peaks(self.dnase_paths[cell_line])
        # get all the chip peaks for this cell line
        chip_peaks = []
        for chip in self.chip_paths[cell_line]:
            chip_peaks = chip_peaks + self.get_peaks(chip)
        # sort the chip peaks by start point then end point
        chip_peaks.sort(key=lambda x: (x[0], x[1]))

        ret = []
        for peak in dnase_peaks:
            ret = ret + self.find_ints(chip_peaks, peak[0], peak[1])
