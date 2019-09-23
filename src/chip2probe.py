import urllib
import os, sys, subprocess, pathlib
import configparser
import datetime
from collections import OrderedDict

from tqdm import tqdm
from collections import defaultdict
from itertools import chain

import pandas as pd

# create virtual environment with python version 2.7 (required by MCAS2)

# we need to import the package folder and libsvm
sys.path.append("probefilter")
sys.path.append("probefilter/libsvm-3.23/python")
# from sitesfinder.imads import iMADS
# from sitesfinder.imadsmodel import iMADSModel
# from sitesfinder.plotcombiner import PlotCombiner
# from sitesfinder.pbmescore import PBMEscore
# from sitesfinder.sequence import Sequence

'''
Summarize
lab-archive -> note the result
information about the data in the plot
'''

# tagsize = 36
# macs_p = 0.01

#bedpath = "/data/gordanlab/vincentius/cooperative_probe/hg19_0005_Ets1.bed"
bedpath = "../imads_files/predictions/hg19_0005_Ets1_filtered.bed"

# Analysis directory
escore_short_path = "../escores/ets1_escores.txt"
escore_map_path = "../escores/index_short_to_long.csv"

# for iMADS, must specify cores and model files
modelcores = ["GGAA",  "GGAT"]
modelpaths = ["../imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
             "../imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model"]
modelwidth = 20 # TODO: confirm if we can get length without manually specify it
imads_cutoff = 0.2128
model_kmers = [1,2,3]

escore_cutoff = 0.4

# ============================

# outdir = "../result/%s" % chipname
init_analysis_file = "sitefiles_list.txt" # this is file obtained after running Rscript

# From https://stackoverflow.com/questions/15644964/python-progress-bar-and-downloads
class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        try: 
            urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)
            return 0
        except:
            return -1


def get_file_info(filename):
    file_info = defaultdict(dict)
    with open(filename, "r") as f:
        next(f)
        for line in f.readlines():
            if line.rstrip() != "" and not line.startswith("#"):
                items = line.strip().split() # split on space and tab
                exp_id, chip_name, rep_tag, file_id, quality, output_type, antibody_id  = \
                items[0], items[1], items[2], items[3], items[4], items[5], items[6]
                if rep_tag.startswith("r"):
                    chip_name = chip_name + "_" + exp_id  
                else:
                    corresponding_chip = items[7]
                    chip_name = chip_name + "_"+ corresponding_chip
                full_tag = rep_tag + "_" + output_type
                file_info[chip_name][full_tag] = file_id
    return file_info


def download_chip(file_info, input_dir, macs_args):

    filtered_file = configparser.ConfigParser()
    unfiltered_file = configparser.ConfigParser()
    outdir = "../result/chipseq"
    timestamp = ""
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for chipname, tag_fileID in file_info.items():
        chip_path = "%s/%s" % (outdir, chipname)
        if not os.path.exists(chip_path):
            os.makedirs(chip_path)
        timestamp += chipname +"\n"
        unfiltered, filtered = OrderedDict(), OrderedDict()
        for tag, file_id in tag_fileID.items():
            fname = file_id + ".bam"
            saveto =  "%s/%s" % (chip_path, fname)
            if not os.path.exists(saveto): 
                chipurl = "https://www.encodeproject.org/files/{0}/@@download/{0}.bam".format(file_id)  
                #print("Downloading %s ..." % (fname))
                status = download_url(chipurl, saveto)  # takes filename?!
                if status == -1:
                    not_download_msg = "%s was NOT downloaded due to an error." % fname
                    timestamp += "\t"+not_download_msg +"\n"
                    print(not_download_msg)
                    continue
                else:
                    now = datetime.datetime.now()
                    timestamp += "\tDownloaded %s on %s at %s \n" % (fname, now.strftime("%Y-%m-%d"), now.strftime("%H:%M:%S"))
                    print("Finished downloading %s" % fname)
            else: # if file already existed
                skip_msg = "%s already existed. Skipped downloading the file." % fname
                timestamp += "\t"+ skip_msg +"\n"
                print(skip_msg)

            # add to config dictionary
            if "unfiltered" in tag:
                unfiltered[tag.split("_")[0]] = fname
            else:
                filtered[tag.split("_")[0]] = fname

            if len(filtered)!=0:
                filtered_file[chipname] = filtered
            if len(unfiltered)!=0:
                unfiltered_file[chipname] = unfiltered

    with open("%s/download_timestamp.txt"% outdir, 'w') as f:
        f.write(timestamp)
    with open("%s/chips_unfiltered.config"% outdir, 'w') as configfile:
        # for call_peaks() to use later
        unfiltered_file.write(configfile)
    with open("%s/chips_filtered.config"%outdir, 'w') as configfile:
        filtered_file.write(configfile)
    
# trying to test fork

def call_peaks(macs_args, config_parser, macs_result_path):
    for chipname in config_parser.sections():
        reps, cntrl = [0]*2, []
        for (tag, filename) in config_parser.items(chipname):
            if tag == "r1": reps[0] = filename
            elif tag == "r2": reps[1] = filename
            elif tag.startswith("c"): cntrl.append(filename)
        all_rep = [[r] for r in reps]
        all_rep.append(reps) # [[r1], [r2], [r1, r2]]
        if 0 in reps: continue
        tagsize, macs_p = macs_args[0], macs_args[1]
        for indx, rep in enumerate(all_rep):
            output_file = ("%s/%s_both"%(macs_result_path, chipname)) if len(rep)==2 else ("%s/%s_r%d"%(macs_result_path, chipname, indx+1))
            rep_str, cntrl_str = ' '.join(rep), ' '.join(cntrl)
            #command = f"macs2 callpeak -t {rep_str} -c {cntrl_str} -n {output_file} -s {tagsize} -p {macs_p} -f BAM -g hs -B"
            command = "macs2 callpeak -t %s -c %s -n %s -s %d -p %f -f BAM -g hs -B" % (rep_str, cntrl_str, output_file, tagsize, macs_p)
            subprocess.call(command, shell=False)
            # print(command)
        

    # analysis_result_path = "%s/analysis_result" % (outdir)
    # if not os.path.exists(analysis_result_path):
    #     os.makedirs(analysis_result_path)
    # print("Running analysis...")
    # pwd = os.path.dirname(os.path.realpath(__file__))
    # pu1_path = "%s/%s%s" % (macs_result_path,chipname,"_r1_treat_pileup.bdg")
    # pu2_path = "%s/%s%s" % (macs_result_path,chipname,"_r2_treat_pileup.bdg")
    # pu_both_path = "%s/%s%s" % (macs_result_path,chipname,"_bothrs_treat_pileup.bdg")
    # nrwp_preidr_path = "%s/%s%s" % (macs_result_path,chipname,"_bothrs_peaks.narrowPeak")
    # nrwp_postidr_path = "%s/%s" % (idr_result_path,"idr_001p_wlist.005i")
    # args_rscript = [pu1_path, pu2_path, pu_both_path, nrwp_preidr_path, nrwp_postidr_path, bedpath, analysis_result_path, chipname]
    # #subprocess.call(["srun","Rscript","R_analysis/main.R",pwd] + args_rscript,shell=False)
    # #subprocess.call(["Rscript","R_analysis/main.R",pwd] + args_rscript,shell=False)



'''
Note: put all inputs in the input dir
Enforce:  
- pipeline.config 

If skipped downloadchip:
- chip_path.config

'''

def main():
    input_dir = "../input/"
    config_pipeline = configparser.ConfigParser()
    config_pipeline.read(input_dir + "pipeline.config") 
    # NOTE: set defaults in pipline, in case users don't follow the template
    
    if(config_pipeline['pipeline'].getboolean('downloadchip')):
        filename = config_pipeline['downloadchip_param']['chip_to_download']
        file = os.path.join(input_dir, filename)
        file_info  = get_file_info(file)
        download_chip(file_info, input_dir)
        print("Finished downloading all files :)")
    else:
        print("Skipping downloadchip!")

    # add try-except to ensure the needed input directories and files are in-place
    if(config_pipeline['pipeline'].getboolean('callpeaks')): # expect to have a 'chipseq' folder in the 'result' folder
        macs_args = (config_pipeline['callpeaks_param']['tagsize'], config_pipeline['callpeaks_param']['macs_p'])
        chip_result_path = "../result/chipseq"
        # filtered & unfiltered - might not be necessary for others
        filtered_file = configparser.ConfigParser()
        unfiltered_file = configparser.ConfigParser()
        filtered_file.read("%s/chips_filtered_test.config"%chip_result_path)
        unfiltered_file.read("%s/chips_unfiltered_test.config"%chip_result_path)
        call_peaks(macs_args, filtered_file, macs_result_path = "../result/called_peaks/filtered")
        #call_peaks(chip_result_path, macs_args, unfiltered_file, macs_result_path = "../result/called_peaks/unfiltered")       
        
        print("Finished calling peaks.")
    else: 
        print("Skipping callpeaks!")

      

if __name__=="__main__":
    main()
    

    

'''

    
    subprocess.call(["srun","macs2.sh",chipdata["r1"],chipdata["r2"],chipdata["c1"],chipdata["c2"],"%s/%s" % (macs_result_path,chipname), str(tagsize)],shell=False)
    #subprocess.call(["./macs2.sh",chipdata["r1"],chipdata["r2"],chipdata["c1"],chipdata["c2"],"%s/%s" % (macs_result_path,chipname), str(tagsize)],shell=False)
    print("Finished running macs, results are saved in %s" % macs_result_path)

    idr_result_path = "%s/idr_result" % (outdir)
    if not os.path.exists(idr_result_path):
        os.makedirs(idr_result_path)
    print("Running idrs...")
    #subprocess.call(["srun","idr.sh","%s/%s" % (macs_result_path,chipname),idr_result_path],shell=False)
    #subprocess.call(["./idr.sh","%s/%s" % (macs_result_path,chipname),idr_result_path],shell=False)


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
    #subprocess.call(["srun","Rscript","R_analysis/main.R",pwd] + args_rscript,shell=False)
    #subprocess.call(["Rscript","R_analysis/main.R",pwd] + args_rscript,shell=False)

    # ============== PLOT AND FILTERING PART ==============

    # First, we can just load the models to avoid having to reload this on every iteration
    models = [iMADSModel(modelpath, modelcore, modelwidth, model_kmers) for modelpath, modelcore in zip(modelpaths, modelcores)]
    imads = iMADS(models, imads_cutoff) # 0.2128 is for the ETS1 cutoff

    escore = PBMEscore(escore_short_path, escore_map_path)

    sitelist_path = "%s/%s" % (analysis_result_path, "sitefiles_list.txt")
    with open(sitelist_path, 'r') as f:
        sitelist = [line.strip() for line in f.readlines()]
    for sitefile in sitelist:
        filename = os.path.basename(os.path.splitext(sitefile)[0])
        print("Making sites plot for %s" % filename)
        seqdf = pd.read_csv(sitefile, sep='\t')

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
            es_pred1 = es_preds[key]
            imads_pred1 =  imads_preds[key]
            bs = BindingSites(imads_pred1,es_pred1)
            if bs.site_count() == 2:
                filtered_sites[key] = bs
        site_list = [{**{"key":site, "sequence":es_preds[site].sequence},**filtered_sites[site].get_sites_dict()} for site in filtered_sites]
        columns = ["key", "site_start_1", "site_start_2", "site_end_1", "site_end_2", "site_pos_1", "site_pos_2", "imads_score_1", "imads_score_2", "sequence"]
        pd.DataFrame(site_list).to_csv(dfpath, index=False, columns=columns, float_format='%.4f')

    print("Done")
'''