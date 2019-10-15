import os, subprocess
def run_idr(macs_result_path, idr_result_path, idr_threshold):
    # loop through dirs
    if not os.path.exists(idr_result_path):
        os.makedirs(idr_result_path)
    for chipname in os.listdir(macs_result_path):
        print(chipname)
        target_dir = macs_result_path+"/"+chipname+"/"+chipname
        r1, r2, both = target_dir+"_r1_peaks.narrowPeak", target_dir+"_r2_peaks.narrowPeak", target_dir+"_both_peaks.narrowPeak"
        out_name, log_name = idr_result_path+"/"+chipname+"_idr.txt", idr_result_path+"/"+chipname+"_idr_log.txt"
        command = "idr --samples %s %s --peak-list %s --idr-threshold %.2f --output-file %s --plot --log-output-file %s" % (r1, r2, both, idr_threshold, out_name, log_name)
        cmd_lst = command.split()
        subprocess.call(cmd_lst, shell=False)

run_idr(macs_result_path = "../result/called_peaks/filtered/to_run_idr", idr_result_path = "../result/idr_output", idr_threshold = 0.05)