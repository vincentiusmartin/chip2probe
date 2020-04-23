pu1_path <- "../../result/ets1_k562/macs_result/ets1_k562_bothreplicates_treat_pileup.bdg"
pu2_path <- "../../result/runx1_k562/macs_result/runx1_k562_bothrs_treat_pileup.bdg"

nrwp1_path <- "../../result/ets1_k562/idr_result/idr_001p_wlist.narrowPeak"
nrwp2_path <- "../../result/runx1_k562/idr_result/idr_001p_wlist.005i"

bed1_path <- "../../resources/imads_files/predictions/hg19_0005_Ets1_filtered.bed"
bed2_path <- "../../resources/imads_files/predictions/hg19_0010_Runx1_filtered.bed"

setwd("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/probe_generator")
source("chip.quality/R/quality_check.R")

peaklen <- 300
allowed_distance <- c(4,24)

genome_sites1 <- read.sites.bed(bed_tf1_path)
genome_sites2 <- read.sites.bed(bed_tf2_path)

pu1 <- read.pileup(pu1_path)
pu2 <- read.pileup(pu2_path)

nrwp1 <- read.narrow.peak(nrwp1_path, peak_length=peaklen, peak_type="summit")
nrwp2 <- read.narrow.peak(nrwp2_path, peak_length=peaklen, peak_type="summit")

agg1 <- calculate.peak.pileup(nrwp1, pu1, logpileup=TRUE) 
agg2 <- calculate.peak.pileup(nrwp2, pu2, logpileup=TRUE) 

# get overlap between peak files
all_sites <- get.heterotypic.sites(agg1, agg2, genome_sites1, genome_sites2, allowed_distance)
all_seqs <- get.probeseq.in.range(all_sites, allowed_distance[1], allowed_distance[2], flank_size = 10)
all_seqs$s1 <- all_seqs$s1 - all_seqs$seqstart
all_seqs$s2 <- all_seqs$s2 - all_seqs$seqstart

write.csv(all_seqs, "test.csv", row.names = FALSE, quote = FALSE)
