#!/bin/bash


# main_dir="chipseq_data/yuning_ets1_k562"
# r1_bam="$main_dir/r1_ENCFF006UXO.bam"
# r2_bam="$main_dir/r2_ENCFF468AKT.bam"
# c1_bam="$main_dir/c1_ENCFF750MIM.bam"
# c2_bam="$main_dir/c2_ENCFF235CSW.bam"
# name2c="ets1_k562"
#tagsize=36
r1_bam=$1
r2_bam=$2
c1_bam=$3
c2_bam=$4
name2c=$5
tagsize=$6


# -s $tagsize
#echo "macs" $r1_bam $r2_bam $c1_bam $c2_bam $name2c $tagsize
macs2 callpeak -t $r1_bam $r2_bam  -c $c1_bam -f BAM -g hs -n ${name2c}_bothrs -B -s $tagsize -p 0.01
macs2 callpeak -t $r1_bam  -c $c1_bam -f BAM -g hs -n ${name2c}_r1 -B -s $tagsize -p 0.01
macs2 callpeak -t $r2_bam  -c $c1_bam -f BAM -g hs -n ${name2c}_r2 -B -s $tagsize -p 0.01
