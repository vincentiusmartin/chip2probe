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


class ProbeFilter:
    """This class generate clean noncustom and custom probes."""

    def __init__(self, **kwargs):
        """Initialize class variables."""
        for k, v in kwargs.items():
            setattr(self, k, v)
        if len(self.colors) == 0:
            self.colors = [('crimson', 'plum'), ('steelblue', 'lightblue')]

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

        self.filter_probes()

        if len(self.customs) > 0:
            self.customize_probes()

    def initialize_objects(self):
        """Initialize escore and model objects."""
        escores = {}
        models = {}
        for protein in self.proteins:
            escores[protein] = PBMEscore(self.escore_long_paths[protein])
            # initialize binding site prediction model used based on method
            if self.method == 'imads':
                # read imads prediction
                imads_models = [iMADSModel(modelpath, modelcore, 20, [1, 2, 3])
                                for modelpath, modelcore in
                                zip(self.modelpaths[protein], self.modelcores[protein])]
                models[protein] = iMADS(imads_models, self.imads_cutoff[protein])
            elif self.method == 'kompas':
                models[protein] = Kompas(protein=protein, threshold=self.mutate_cutoff,
                                         kmer_align_path=self.kmer_align_paths[protein])
            # raise exception if model is not supported
            else:
                raise Exception("Supported methods are imads and kompas. Got {}".format(method))
        self.escores = escores
        self.models = models

    def filter_probes(self):
        """Filter and clean up probes so each sequence has only 2 significant sites."""
        sitepath = "%ssites_all.tsv" % self.analysis_path
        seqdf = pd.read_csv(sitepath, sep='\t')

        # get all combinations of number of sites in peak and peak length as a list of tuples
        spcomb = sorted(list(set(zip(seqdf["sites_in_peak"], seqdf["peaklen"]))), key=lambda x: (x[0],x[1]))

        # for each combination, get filtered probes
        spcomb = [(3, 200)]

        filtered_probes = coopfilter.get_filtered_probes(seqdf=seqdf,
                                                         escores=self.escores,
                                                         models=self.models,
                                                         mutate_cutoff=self.mutate_cutoff,
                                                         mutate_gap=self.mutate_gap,
                                                         egaps=self.egaps,
                                                         thresholds=self.thresholds,
                                                         proteins=self.proteins,
                                                         colors=self.colors,
                                                         generate_plots=self.generate_plots,
                                                         spcomb=spcomb,
                                                         analysis_path=self.analysis_path,
                                                         mode="noncustom",
                                                         predict_flanks=True,
                                                         key_colname="key",
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
            print("Saving result in %smutated_probes.tsv" % self.analysis_path)
            df_out.to_csv("%smutated_probes.tsv" % (self.analysis_path), sep="\t", index=False)
        else:
            print("No mutation rows found in %s" % sitepath)

    def customize_probes(self):
        """Customize clean probes."""
        infile = "%smutated_probes.tsv" % self.analysis_path
        df = pd.read_csv(infile, sep="\t")
        df = df[:10]
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
