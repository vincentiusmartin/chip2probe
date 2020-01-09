#!/bin/bash

#name2c="macs_result/yuning_ets1_k562_2/ets1_k562"
chipname=$1
outpath=$2

#echo "idr" $name2c $outpath
idr -s ${chipname}_r1_peaks.narrowPeak ${chipname}_r2_peaks.narrowPeak -p ${chipname}_both_peaks.narrowPeak --input-file-type narrowPeak --plot -i 0.05 -o $outpath/idr_001p_wlist.005i
