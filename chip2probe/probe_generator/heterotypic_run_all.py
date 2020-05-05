import os
import argparse
import urllib.request
import subprocess

from tqdm import tqdm

from heterotypic_dictinput_example import inlist

# todo: fix the merge similar function with homotypic

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)

def handle_encode_url(indict, outdirs):
    result_dir = {}
    # iterate over 2 transcription factors
    for i in range(len(outdirs)):
        tfid =  "tf%d" % (i+1)
        tfname = indict["%s_name" % tfid]
        chipdata_path = "%s/chipseq_data" % (outdirs[i])
        if not os.path.exists(chipdata_path):
            os.makedirs(chipdata_path)
        chip_info = "ChIP-seq data for %s:\n" % indict["%s_name" % tfid]

        chip_list = ["%s_r1"%tfid,"%s_r2"%tfid,"%s_c1"%tfid,"%s_c2"%tfid]
        #get available replicates and controls
        files_dl = [ch for ch in chip_list if ch in indict.keys()]
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
        macs_result_path = "%s/macs_result" % (outdirs[i])
        if not os.path.exists(macs_result_path):
            os.makedirs(macs_result_path)
        print("Running macs...")
        chip_files = [chipdata[x] if x in indict.keys() else "-" for x in chip_list]
        macs_cmdlist = ["./macs2.sh"] + chip_files + ["%s/%s" % (macs_result_path,tfname), str(tagsize)]
        print(" ".join(macs_cmdlist))
        #subprocess.call(macs_cmdlist,shell=False)
        print("Finished running macs, results are saved in %s" % macs_result_path)

        idr_result_path = "%s/idr_result" % (outdirs[i])
        if not os.path.exists(idr_result_path):
            os.makedirs(idr_result_path)
        print("Running idr...")
        idr_cmdlist = ["./idr.sh","%s/%s" % (macs_result_path,tfname),idr_result_path]
        print(" ".join(idr_cmdlist))
        #subprocess.call(idr_cmdlist,shell=False)
        print("Finished running idr, results are saved in %s" % idr_result_path)

        result_dir[tfid] = {}
        result_dir[tfid]["pu1_path"] = "%s/%s%s" % (macs_result_path,tfname,"_r1_treat_pileup.bdg")
        result_dir[tfid]["pu2_path"] = "%s/%s%s" % (macs_result_path,tfname,"_r2_treat_pileup.bdg")
        result_dir[tfid]["pu_both_path"] = "%s/%s%s" % (macs_result_path,tfname,"_bothreplicates_treat_pileup.bdg")
        result_dir[tfid]["nrwp_preidr_path"] = "%s/%s%s" % (macs_result_path,tfname,"_bothreplicates_peaks.narrowPeak")
        result_dir[tfid]["nrwp_postidr_path"] = "%s/%s" % (idr_result_path,"idr_001p_wlist.narrowPeak")

    return result_dir


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Run the pipeline.')
    parser.add_argument('index', type=int, help="Input index from dictinput.py")
    parser.add_argument('-o','--outdir', type=str, help="directory to output the results", default=".")
    args = parser.parse_args()

    findex = args.index
    # we still use the individual directory of each tf, but create another one that
    # has both
    outdirs = [
        "%s/%s" % (args.outdir, inlist[findex]["tf1_name"]),
        "%s/%s" % (args.outdir, inlist[findex]["tf2_name"])
    ]
    heterodir = "%s/%s_%s" % (args.outdir, inlist[findex]["tf1_name"], inlist[findex]["tf2_name"])

    for outdir in outdirs + [heterodir]:
        if not os.path.exists(outdir):
            print("Make directory %s" % outdir)
            os.makedirs(outdir)

    # if it is from encode, we need to download it and get the required files
    if inlist[findex]["type"] == "encode":
        tfdir = handle_encode_url(inlist[findex], outdirs)
        args_rscript = [
            tfdir["tf1"]["pu_both_path"],
            tfdir["tf1"]["nrwp_postidr_path"],
            tfdir["tf2"]["pu_both_path"],
            tfdir["tf2"]["nrwp_postidr_path"]
        ]
    elif inlist[findex]["type"] == "cistrome":
        args_rscript = [
            "-",
            inlist[findex]["tf1_path"],
            "-",
            inlist[findex]["tf2_path"]
        ]
    else:
        raise Exception("Type is not supported")
    # skip cistrome since the data doesn't need processing

    analysis_result_path = "%s/analysis_result" % heterodir
    if not os.path.exists(analysis_result_path):
        os.makedirs(analysis_result_path)
    print("Getting heterotypic sites...")
    pwd = os.path.dirname(os.path.realpath(__file__))
    qc_cmdlist = ["Rscript","chip.quality/main_heterotypic.R",pwd] + args_rscript + [analysis_result_path]
    print(" ".join(qc_cmdlist))
    subprocess.call(qc_cmdlist,shell=False)

    # run probefilter here
