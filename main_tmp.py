import urllib.request
import os
import subprocess
import pandas as pd

from tqdm import tqdm

import sys
sys.path.append("probefilter")
sys.path.append("probefilter/libsvm-3.23/python")
from sitesfinder.imads import iMADS
from sitesfinder.imadsmodel import iMADSModel
from sitesfinder.plotcombiner import PlotCombiner
from sitesfinder.pbmescore import PBMEscore
from sitesfinder.sequence import Sequence
from sitesfinder.prediction.basepred import BasePrediction
from cooperative import coopfilter

'''
Summarize
lab-archive -> note the result
information about the data in the plot
'''

chipname = "ets1_HepG2"
chipurls = {
    "r1":"https://www.encodeproject.org/files/ENCFF477EHC/@@download/ENCFF477EHC.bam",
    "r2":"https://www.encodeproject.org/files/ENCFF371ZBY/@@download/ENCFF371ZBY.bam",
    "c1":"https://www.encodeproject.org/files/ENCFF963CVB/@@download/ENCFF963CVB.bam",
    #"c2":"https://www.encodeproject.org/files/ENCFF956QBZ/@@download/ENCFF956QBZ.bam",
}
tagsize = 36

#bedpath = "/data/gordanlab/vincentius/cooperative_probe/hg19_0005_Ets1.bed"
bedpath = "hg19_0005_Ets1.bed"

# Analysis directory
escore_short_path = "../escores/ets1_escores.txt"
escore_map_path = "../escores/index_short_to_long.csv"

# for iMADS, must specify cores and model files
modelcores = ["GGAA",  "GGAT"]
modelpaths = ["../imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
             "../imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"]
modelwidth = 20 # TODO: confirm if we can get length without manually specifying it
imads_cutoff = 0.2128
model_kmers = [1,2,3]

escore_cutoff = 0.4

# ============================

outdir = "../result/%s" % chipname

# From https://stackoverflow.com/questions/15644964/python-progress-bar-and-downloads
class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)

if __name__=="__main__":
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    chipdata_path = "%s/chipseq_data" % (outdir)
    if not os.path.exists(chipdata_path):
        os.makedirs(chipdata_path)
    chipdata = {}
    chip_info = "ChIP-seq data for %s:\n" % chipname
    # ===== Download ChIP-seq data =====
    for key in chipurls:
        fname = os.path.basename(chipurls[key])
        saveto = os.path.join(chipdata_path, fname)
        chipdata[key] = saveto
        chip_info += "%s: %s\n" % (key,fname)
        print("Downloading %s to %s:" % (key,saveto))
        #download_url(chipurls[key], saveto)
    with open("%s/chipinfo.txt" % (outdir), 'w') as f:
        f.write(chip_info)

    macs_result_path = "%s/macs_result" % (outdir)
    if not os.path.exists(macs_result_path):
        os.makedirs(macs_result_path)
    print("Running macs...")
    #subprocess.call(["srun","macs2.sh",chipdata["r1"],chipdata["r2"],chipdata["c1"],"%s/%s" % (macs_result_path,chipname), str(tagsize)],shell=False)
    print("Finished running macs, results are saved in %s" % macs_result_path)

    idr_result_path = "%s/idr_result" % (outdir)
    if not os.path.exists(idr_result_path):
        os.makedirs(idr_result_path)
    print("Running idrs...")
    #subprocess.call(["srun","idr.sh","%s/%s" % (macs_result_path,chipname),idr_result_path],shell=False)

    analysis_result_path = "%s/analysis_result" % (outdir)
    if not os.path.exists(analysis_result_path):
        os.makedirs(analysis_result_path)
    print("Running analysis...")
    pwd = os.path.dirname(os.path.realpath(__file__))
    pu1_path = "%s/%s%s" % (macs_result_path,chipname,"_r1_treat_pileup.bdg")
    pu2_path = "%s/%s%s" % (macs_result_path,chipname,"_r2_treat_pileup.bdg")
    pu_both_path = "%s/%s%s" % (macs_result_path,chipname,"_bothrs_treat_pileup.bdg")
    nrwp_preidr_path = "%s/%s%s" % (macs_result_path,chipname,"_bothrs_peaks.narrowPeak")
    nrwp_postidr_path = "%s/%s" % (idr_result_path,"idr_001p_wlist.005i")
    args_rscript = [pu1_path, pu2_path, pu_both_path, nrwp_preidr_path, nrwp_postidr_path, bedpath, analysis_result_path, chipname]
    #print(["R_analysis/main.R",pwd] + args_rscript)
    #subprocess.call(["srun","Rscript","R_analysis/main.R",pwd] + args_rscript,shell=False)
    #subprocess.call(["Rscript","R_analysis/main.R",pwd + "/../"] + args_rscript,shell=False)

    # ============== PLOT AND FILTERING PART ==============

    # First, we can just load the models to avoid having to reload this on every iteration
    models = [iMADSModel(modelpath, modelcore, modelwidth, model_kmers) for modelpath, modelcore in zip(modelpaths, modelcores)]
    imads = iMADS(models, imads_cutoff) # 0.2128 is for the ETS1 cutoff

    escore = PBMEscore(escore_short_path, escore_map_path)

    sitelist_path = "%s/%s" % (analysis_result_path, "sitefiles_list.txt")
    with open(sitelist_path, 'r') as f:
        sitelist = [line.strip() for line in f.readlines()]
    for sitepath in sitelist:
        print("dsa",sitepath)
        filename = os.path.basename(os.path.splitext(sitepath)[0])
        print("Making sites plot for %s" % filename)
        seqdf = pd.read_csv(sitepath, sep='\t')

        # Make Escore object
        es_preds = escore.predict_sequences(seqdf)
        eplots = escore.plot(es_preds)
        # Make iMADS plot
        imads_preds = imads.predict_sequences(seqdf)
        imadsplots = imads.plot(imads_preds)
        plots = [imadsplots, eplots]

        pc = PlotCombiner() # can do this just once but not a big deal
        plotpath = "%s/sitesplot_%s.pdf" % (analysis_result_path, filename)
        pc.plot_seq_combine(plots, filepath=plotpath)

        filtered_sites = {}
        print("Site filtering...")
        for key in es_preds:
            bs = Sequence(es_preds[key],imads_preds[key])
            if bs.site_count() == 2:
                filtered_sites[key] = bs
        #site_list = [{**{"key":site, "sequence":es_preds[site].sequence},**filtered_sites[site].get_sites_dict()} for site in filtered_sites]
        #columns = ["key", "site_start_1", "site_start_2", "site_end_1", "site_end_2", "site_pos_1", "site_pos_2", "imads_score_1", "imads_score_2", "sequence"]
        #pd.DataFrame(site_list).to_csv("%s/sitelist_%s.pdf" % (analysis_result_path), index=False, columns=columns, float_format='%.4f')

        seqdict = {}
        funcdict = {}
        filtered_probes = []
        # TODO: tmr look at 110,271
        for key in filtered_sites:
        #for key in ["sequence11"]:
            # Visualization part
            seqdict["%s-wt" % key] = filtered_sites[key].sequence
            for idx,mut in enumerate([[0],[1],[0,1]]):
                mutseq = filtered_sites[key].abolish_sites(mut,escore)
                seqdict["%s-m%d" % (key,idx + 1)] = mutseq.sequence
                funcdict["%s-m%d" % (key,idx + 1)] = mutseq.plot_functions
            if coopfilter.filter_coopseq(seqdict["%s-wt"%key], seqdict["%s-m1"%key],
                                         seqdict["%s-m2"%key], seqdict["%s-m3"%key],
                                         filtered_sites[key].get_sites_dict(), escore):
                filtered_probes.append({"key":key, "wt":seqdict["%s-wt"%key], "m1":seqdict["%s-m1"%key],
                                        "m2":seqdict["%s-m2"%key], "m3":seqdict["%s-m3"%key]})

        pp = escore.plot(escore.predict_sequences(seqdict),additional_functions=funcdict)
        pc.plot_seq_combine([pp], filepath="%s/plot_mut_%s.pdf" % (analysis_result_path,filename))

        pd.DataFrame(filtered_probes).to_csv("%s/mutated_probes_%s.tsv" % (analysis_result_path,filename),sep="\t",index=False,columns=["key","wt","m1","m2","m3"])

    #print(fname,header)
