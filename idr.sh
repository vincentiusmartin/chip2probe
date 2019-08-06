#!/bin/bash

#name2c="macs_result/yuning_ets1_k562_2/ets1_k562"
name2c=$1
outpath=$2

#echo "idr" $name2c $outpath
idr -s ${name2c}_r1_peaks.narrowPeak ${name2c}_r2_peaks.narrowPeak -p ${name2c}_bothrs_peaks.narrowPeak --input-file-type narrowPeak --plot -i 0.05 -o $outpath/idr_001p_wlist.005i
