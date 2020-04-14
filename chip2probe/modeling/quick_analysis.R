library(dplyr)

setwd("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/modeling")

df <- read.csv("/Users/vincentiusmartin/Research/chip2gcPBM/probedata/191004_coop-PBM_Ets1_v1_1st/training_data/training.tsv", sep ="\t")

grouped <- df %>% group_by(label)
