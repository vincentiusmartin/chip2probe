
# these paths need to be handled better
import sys
sys.path.append('libsvm-3.23/python')
sys.path.append('..')
sys.path.append("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe")
from util import bio

from sitespredict.pbmescore import PBMEscore
from sitespredict.imadsmodel import iMADSModel
from sitespredict.sitesplotter import SitesPlotter
from sitespredict.sequence import Sequence
from sitespredict.imads import iMADS

from customseq import customseq

import pandas as pd
import pickle

# TODO TODO
def coop_probes_to_df(result_path, pattern):
    re_pattern = re.compile(pattern)
    agg = []
    for path, dirs, filenames in os.walk(result_path):
        for fname in filter(lambda name:re_pattern.match(name),filenames):
            tf_name = path.split("/")[-2] # should be the second from the last
            fsplit = os.path.splitext(fname)[0].split("_")
            span = fsplit[-1][len("span"):]
            nsites = fsplit[-2][1:]
            filepath = os.path.join(path, fname)

            probe_df = pd.read_csv(filepath, sep="\t")
            probe_df["tf"] = tf_name
            probe_df["span"] = span
            probe_df["nsites"] = nsites
            agg.extend(probe_df.to_dict('records'))
    df = pd.DataFrame(agg, columns=["wt","m1","m2","m3","ecutoff","nsites","egapthres","span","tf","flank_left","flank_right","key"], index=None)
    df = df.rename(columns={'key': 'seqlabel'})
    return df

def seqdf2mutfile(df, escore, imads,key_colname="key"):
    """
    flank_left and flank_right, both are dictionaries
    """

    print("Start making dist sites file...")

    flank_left = bio.get_seqdict(df,"flank_left",keycolname="key")
    flank_right = bio.get_seqdict(df,"flank_right",keycolname="key")

    # Make Escore object
    es_preds = escore.predict_sequences(df,key_colname="key")
    # Make iMADS plot
    imads_preds = imads.predict_sequences(df,key_colname="key")

    filtered_sites = {}
    flanks = {}
    print("Site filtering...")
    for key in es_preds:
        bs = Sequence(es_preds[key],imads_preds[key],escore_cutoff=0.2) # 0.2 for distance only
        if bs.site_count() == 2:
            filtered_sites[key] = bs

    seqdict = {}
    funcdict = {}
    filtered_probes = []
    for key in filtered_sites:
    #for key in ["sequence11"]:
        # Visualization part
        curseqdict = {}
        curfuncdict = {}
        curseqdict["%s-wt" % key] = filtered_sites[key].sequence
        for idx,mut in enumerate([[0],[1],[0,1]]):
            mutseq = filtered_sites[key].abolish_sites(mut,escore)
            curseqdict["%s-m%d" % (key,idx + 1)] = mutseq.sequence
            curfuncdict["%s-m%d" % (key,idx + 1)] = mutseq.plot_functions
        if len(set(curseqdict.values())) == 4:
            seqdict.update(curseqdict)
            funcdict.update(curfuncdict)
            bsites_dict = filtered_sites[key].get_sites_dict()
            filtered_probes.append({"key":key,
                                "wt":seqdict["%s-wt"%key],
                                "m1":seqdict["%s-m1"%key],
                                "m2":seqdict["%s-m2"%key],
                                "m3":seqdict["%s-m3"%key],
                                "flank_left":flank_left[key],
                                "flank_right":flank_right[key],
                                "core1_start":bsites_dict["core_start_1"],
                                "core1_end":bsites_dict["core_end_1"],
                                "site1_pref":bsites_dict["imads_score_1"],
                                "core2_start":bsites_dict["core_start_2"],
                                "core2_end":bsites_dict["core_end_2"],
                                "site2_pref":bsites_dict["imads_score_2"],
                                "ecutoff": 0.2,
                                "egapthres":"-",
                                "distance":filtered_sites[key].get_sites_dist(),
                                "sites_in_peak": "-",
                                "peak_length": "-",
                                "coordinate": "-"
                            })
        else:
            print("Couldn't mutate sequences in %s" % key)

    #sp = SitesPlotter()
    #pp = escore.plot(escore.predict_sequences(seqdict),additional_functions=funcdict)
    #pc.plot_seq_combine([pp], filepath="plot_distsites.pdf")

    # probably should check here if filtered_probes is empty
    fp_df = pd.DataFrame(filtered_probes)
    req_cols = ["key", "wt", "m1", "m2", "m3", "flank_left", "flank_right", "core1_start", "core1_end", "site1_pref",
                "core2_start", "core2_end", "site2_pref", "ecutoff", "egapthres", "distance"]
    fp_df.to_csv("mutated_probes_dist_and_weak.tsv",columns=req_cols,index=False)

if __name__ == "__main__":
    ncustom = 750 # how many custom seq we want?

    # this is where we parse the directory and get all files
    infile = "/Users/vincentiusmartin/Research/chip2gcPBM/result/ets1_k562/analysis_result/mutated_probes.tsv"
    df = pd.read_csv(infile,sep="\t")
    tf = "ets1_k562"

    # read escore file
    escore_long_path = "/Users/vincentiusmartin/Research/chip2gcPBM/resources/escores/Ets1_8mers_11111111.txt"
    escore = PBMEscore(escore_long_path)

    modelcores = ["GGAA",  "GGAT"]
    modelpaths = ["/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
    "/Users/vincentiusmartin/Research/chip2gcPBM/resources/imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"]
    # read imads prediction
    imads_models = [iMADSModel(modelpath, modelcore, 20, [1,2,3]) for modelpath, modelcore in zip(modelpaths, modelcores)]
    imads = iMADS(imads_models, 0.2128) # 0.2128 is for the ETS1 cutoff

    df_wseq = pd.Series(df["flank_left"].values + df["wt"].values + df["flank_right"].values,index=tf + "_" + df["key"]).to_dict()

    weak_custom = []

    # here, assume sequence have 2 sites
    for key in df_wseq:
        print("Working on %s" % key)
        seq = df_wseq[key]
        print("  Making weak site for site1")
        w1 = customseq.make_weak_site(seq, escore, imads, which_site = 0, mutable_flank_length = 2, left_outer = True) # left
        if seq == w1:
            print("could not find weak sites for the first site")
            continue
        print("  Making weak site for site2")
        w2 = customseq.make_weak_site(seq, escore, imads, which_site = 1, mutable_flank_length = 2, left_outer = False) # right
        if seq == w2:
            print("could not find weak sites for the second site")
            continue
        print("  Making weak site for site1 and site2")
        w3 = customseq.make_weak_site(w1, escore, imads, which_site = 1, mutable_flank_length = 2, left_outer = False) # both
        if w1 == w3 or seq == w3:
            print("could not find weak sites for both sites")
            continue
        weak_custom.append({"key":"%s_weak_s1" % key,"flank_left":w1[0:10],"sequence":w1[10:46],"flank_right":w1[46:56]})
        weak_custom.append({"key":"%s_weak_s2" % key,"flank_left":w2[0:10],"sequence":w2[10:46],"flank_right":w2[46:56]})
        weak_custom.append({"key":"%s_weak_s1_2" % key,"flank_left":w3[0:10],"sequence":w3[10:46],"flank_right":w3[46:56]})
        # if we have enough then break
        if len(weak_custom) >= ncustom * 3:
            break

    # Second, play with the position of the binding sites
    dist_custom = []
    for key in df_wseq:
        print("Working on %s" % key)
        seq = df_wseq[key]
        ims_pred = imads.predict_sequence(seq)
        if len(ims_pred) != 2:
            print("Number of sites is not equal to 2 in %s" % key)
            continue
        # minus 10 because of the flank length
        site1_start = ims_pred[0]["core_start"] - 10
        site2_start = ims_pred[1]["core_start"] - 10
        if site2_start - site1_start < 20:
            print("Distance is less than 20 for %s" % key)
            continue
        site1_end = site1_start + ims_pred[0]["core_width"]
        site2_end = site2_start + ims_pred[1]["core_width"]

        flank_left = df_wseq[key][0:10]
        seq = df_wseq[key][10:46]
        flank_right = df_wseq[key][46:56]
        step = 0 # how much apart should we move the sites
        delta = site2_start - site1_start
        while delta > 5:
            step += 1
            s1 = customseq.move_bsite(seq, site1_start, site1_end, flank_left, step_to_right = step)
            s2 = customseq.move_bsite(s1, site2_start, site2_end, flank_right, step_to_right = -step)
            delta = site2_start - site1_start - step * 2
            dist_custom.append({"key":"%s_dist_%d" % (key,delta), "flank_left":flank_left, "sequence":s2, "flank_right":flank_right})

    #pickle.dump(weak_custom + dist_custom, open("test.pickle","wb"))
    #df_comb = pd.DataFrame(pickle.load(open("test.pickle","rb")))

    df_comb = pd.DataFrame(weak_custom + dist_custom)
    seqdf2mutfile(df_comb, escore, imads)
