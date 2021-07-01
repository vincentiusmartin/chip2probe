args <- commandArgs(trailingOnly=TRUE)
print(args)
setwd(args[1])
pu_tf1_path <- args[2]
nrwp_tf1_path <- args[3]
pu_tf2_path <- args[4]
nrwp_tf2_path <- args[5]
outpath <- args[6]

setwd("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/chip2probe/probe_generator")
pu_tf1_path <- "/Users/vincentiusmartin/Research/chip2gcPBM/result/ets1_k562/macs_result/ets1_k562_bothreplicates_treat_pileup.bdg"
nrwp_tf1_path <- "/Users/vincentiusmartin/Research/chip2gcPBM/result/ets1_k562/idr_result/idr_001p_wlist.narrowPeak"
pu_tf2_path <- "/Users/vincentiusmartin/Research/chip2gcPBM/result/runx1_k562/macs_result/runx1_k562_bothreplicates_treat_pileup.bdg"
nrwp_tf2_path <- "/Users/vincentiusmartin/Research/chip2gcPBM/result/runx1_k562/idr_result/idr_001p_wlist.narrowPeak"

source("chip.quality/R/quality_check.R")

peaklen <- 300
allowed_distance <- c(4,24)

# TODO: put peak pileup in the result
config <- config::get()
genome_sites1 <- get.genome.sites(config, config$tf[1], nrwp_tf1_path, genomever=config$genomever)
genome_sites2 <- get.genome.sites(config, config$tf[2], nrwp_tf2_path, genomever=config$genomever)

nrwp1 <- read.narrow.peak(nrwp_tf1_path, peak_length=peaklen, peak_type="summit")
if (pu_tf1_path != "-"){
  pu1 <- read.pileup(pu_tf1_path)
  agg1 <- calculate.peak.pileup(nrwp1, pu1, logpileup=TRUE) 
}else{ # if pileup is not available, just use the narrow peak file
  agg1 <- nrwp1 %>% dplyr::rename(
    peak_start = start,
    peak_end = end
  )
  setDT(agg1, key=c("chr", "peak_start", "peak_end"))
}

nrwp2 <- read.narrow.peak(nrwp_tf2_path, peak_length=peaklen, peak_type="summit")
if (pu_tf2_path != "-"){
  pu2 <- read.pileup(pu_tf2_path)
  agg2 <- calculate.peak.pileup(nrwp2, pu2, logpileup=TRUE) 
}else{
  agg2 <- nrwp2 %>% dplyr::rename(
    peak_start = start,
    peak_end = end
  )
  setDT(agg2, key=c("chr", "peak_start", "peak_end"))
}


# get overlap between peak files
all_sites <- get.heterotypic.sites(agg1, agg2, genome_sites1, genome_sites2, allowed_distance)
all_seqs <- get.probeseq.in.range(all_sites, allowed_distance[1], allowed_distance[2], flank_size = 10)
all_seqs$s1 <- all_seqs$mid_1 - all_seqs$seqstart - 1
all_seqs$s2 <- all_seqs$mid_2 - all_seqs$seqstart - 1

all_seqs$key <- sprintf("sequence%s",seq.int(nrow(all_seqs)))
tbl_path <- paste(outpath,"/sites_all.tsv",sep='')
cat("Write result to",tbl_path)
write.table(all_seqs,file=tbl_path,sep="\t",row.names = FALSE, quote = FALSE)
