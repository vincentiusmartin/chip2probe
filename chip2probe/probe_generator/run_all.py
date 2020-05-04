import sys
import urllib.request
import os
import subprocess
import pandas as pd
import itertools
import argparse
import yaml

from tqdm import tqdm

from dictinput import inlist

#read the config file here
with open("config.yml", "r") as ymlfile:
    conf = yaml.load(ymlfile, Loader=yaml.FullLoader)['default']

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)

def handle_encode_url(indict, outdir):
    chipdata_path = "%s/chipseq_data" % (outdir)
    if not os.path.exists(chipdata_path):
        os.makedirs(chipdata_path)
    chip_info = "ChIP-seq data for %s:\n" % indict["name"]

    chip_list = ["r1","r2","c1","c2"]
    #get available replicates and controls
    files_dl = [x for x in chip_list if x in indict.keys()]
    # ===== Download ChIP-seq data =====
    chipdata = {}
    for key in files_dl:
        fname = os.path.basename(indict[key])
        saveto = os.path.join(chipdata_path, fname)
        chipdata[key] = saveto
        chip_info += "%s: %s\n" % (key,fname)
        print("Downloading %s to %s:" % (key,saveto))
        #download_url(indict[key], saveto)

    # ===== Running MACS2 =====
    tagsize = 36 # default tagsize
    macs_result_path = "%s/macs_result" % (outdir)
    if not os.path.exists(macs_result_path):
        os.makedirs(macs_result_path)
    print("Running macs...")
    chip_files = [chipdata[x] if x in indict.keys() else "-" for x in chip_list]
    macs_cmdlist = ["./macs2.sh"] + chip_files + ["%s/%s" % (macs_result_path,indict["name"]), str(tagsize)]
    print(" ".join(macs_cmdlist))
    #subprocess.call(macs_cmdlist,shell=False)
    print("Finished running macs, results are saved in %s" % macs_result_path)

    idr_result_path = "%s/idr_result" % (outdir)
    if not os.path.exists(idr_result_path):
        os.makedirs(idr_result_path)
    print("Running idr...")
    idr_cmdlist = ["./idr.sh","%s/%s" % (macs_result_path,indict["name"]),idr_result_path]
    print(" ".join(idr_cmdlist))
    #subprocess.call(idr_cmdlist,shell=False)
    print("Finished running idr, results are saved in %s" % idr_result_path)

    analysis_result_path = "%s/analysis_result" % (outdir)
    if not os.path.exists(analysis_result_path):
        os.makedirs(analysis_result_path)
    print("Running analysis...")
    pwd = os.path.dirname(os.path.realpath(__file__))
    pu1_path = "%s/%s%s" % (macs_result_path,indict["name"],"_r1_treat_pileup.bdg")
    pu2_path = "%s/%s%s" % (macs_result_path,indict["name"],"_r2_treat_pileup.bdg")
    pu_both_path = "%s/%s%s" % (macs_result_path,indict["name"],"_bothreplicates_treat_pileup.bdg")
    nrwp_preidr_path = "%s/%s%s" % (macs_result_path,indict["name"],"_bothreplicates_peaks.narrowPeak")
    nrwp_postidr_path = "%s/%s" % (idr_result_path,"idr_001p_wlist.narrowPeak")
    args_rscript = [pu1_path, pu2_path, pu_both_path, nrwp_preidr_path, nrwp_postidr_path, analysis_result_path, indict["name"]]
    qc_cmdlist = ["Rscript","chip.quality/main.R",pwd] + args_rscript
    print(" ".join(qc_cmdlist))
    subprocess.call(qc_cmdlist,shell=False)

def handle_cistrome(indict, outdir):
    analysis_result_path = "%s/analysis_result" % (outdir)
    if not os.path.exists(analysis_result_path):
        os.makedirs(analysis_result_path)
    print("Running analysis...")
    pwd = os.path.dirname(os.path.realpath(__file__))
    args_rscript = [indict["path"], analysis_result_path, indict["name"]]
    qc_cmdlist = ["Rscript","chip.quality/main_cistrome.R",pwd] + args_rscript
    print(" ".join(qc_cmdlist))
    subprocess.call(qc_cmdlist,shell=False)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Run the pipeline.')
    parser.add_argument('index', type=int, help="Input index from dictinput.py")
    parser.add_argument('-o','--outdir', type=str, help="directory to output the results", default=".")
    args = parser.parse_args()

    findex = args.index
    outdir = "%s/%s" % (args.outdir,inlist[findex]["name"])


    if not os.path.exists(outdir):
        print("Make directory %s" % outdir)
        os.makedirs(outdir)

    # if it is from encode, we need to download it and get the required files
    if inlist[findex]["type"] == "encode":
        handle_encode_url(inlist[findex], outdir)
    elif inlist[findex]["type"] == "cistrome":
        handle_cistrome(inlist[findex], outdir)
    else:
        raise Exception("Type is not supported")


    """
    # We need sites_all.tsv from the previous part, in here, we should have the
    # file already.
    # TODO: we need to avoid hardcode appending path like this
    pwd = os.path.dirname(os.path.realpath(__file__))
    pfilter_cmd = ["python3", "probefilter/main.py", pwd, "%s/analysis_result" % outdir]
    print(" ".join(pfilter_cmd))
    subprocess.call(pfilter_cmd,shell=False)
    #ProbeFilter(**input)._run_all()
    """
