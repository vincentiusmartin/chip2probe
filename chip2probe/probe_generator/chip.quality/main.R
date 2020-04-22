# args <- commandArgs(trailingOnly=TRUE)
# print(args)
# setwd(args[1])
# pu1_path <- args[2]
# pu2_path <- args[3]
# pu_both_path <- args[4]
# nrwp_preidr_path <- args[5]
# nrwp_postidr_path <- args[6]
# bed_path <- args[7]
# outpath <- args[8]
# chip_name <- args[9]

setwd("/Users/vincentiusmartin/Research/chip2gcPBM")
pu1_path <- "result/ets1_k562/macs_result/ets1_k562_r1_treat_pileup.bdg"
pu2_path <- "result/ets1_k562/macs_result/ets1_k562_r2_treat_pileup.bdg"
pu_both_path <- "result/ets1_k562/macs_result/ets1_k562_bothreplicates_treat_pileup.bdg"
nrwp_preidr_path <- "result/ets1_k562/macs_result/ets1_k562_bothreplicates_peaks.narrowPeak"
nrwp_postidr_path <- "result/ets1_k562/idr_result/idr_001p_wlist.narrowPeak"
bed_path <- "resources/imads_files/predictions/hg19_0005_Ets1_filtered.bed" # change var to imads_path
outpath <- "result/ets1_k562/analysis_result"
chip_name <- "ets1_k562"
tfname <- "ets1"

source("chip2probe/chip2probe/probe_generator/chip.quality/R/quality_check.R")

probe_size <- 36
probeseq_flank <- 10 # n to the left and n to the right
peaklen <- c(100,200,300)
count_sites_per_peak <- c(2,3,4)
min_site_dist <- 1
max_site_dist <- 24

# ===== END OF CONFIGURATION ===== 

cat("Reading binding sites prediction...\n")
genome_sites <- read.sites.bed(bed_path)
cat("Reading pileup files...\n")
pu1 <- read.pileup(pu1_path)
pu2 <- read.pileup(pu2_path)
pu_both <- read.pileup(pu_both_path)

# first, get general information about the peak
nrwp_preidr <- read.narrow.peak(nrwp_preidr_path, peak_type="all", show_peaklen=TRUE)
nrwp_postidr <- read.narrow.peak(nrwp_postidr_path, peak_type="all", show_peaklen=TRUE)
write.table(nrwp_preidr,file=paste(outpath,"/narrowpeak_preidr.tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
write.table(nrwp_postidr,file=paste(outpath,"/narrowpeak_postidr.tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
plot.peaklen.dist(list(nrwp_preidr$peaklen, nrwp_postidr$peaklen), c("before_idr", "after_idr"), paste(outpath,"/peaklen_dist.pdf",sep=''), chip_name)

# ---

datalist <- list()
idx <- 0
for (len in peaklen){
  cat("Working on analysis for peaklen",len,"...\n")
  nrwp <- read.narrow.peak(nrwp_postidr_path, peak_type = "summit", peak_length=len)
  peaklenlabel <- paste(chip_name,", peaklen ",len, sep='')
  
  cat("  Generating correlation plot...\n")
  agg1 <- calculate.peak.pileup(nrwp, pu1, logpileup=FALSE)
  agg2 <- calculate.peak.pileup(nrwp, pu2, logpileup=FALSE)
  merged <- merge(x = agg1, y = agg2, by = c("chr", "peak_start", "peak_end"))
  ###write.table(merged,file=paste(outpath,"/pileup_scores_span",span,".tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
  corr_plot_path <- paste(outpath,"/corr_plot_",len,".pdf",sep='')
  
  sites_peak <- foverlaps(genome_sites,nrwp, nomatch=NULL)[, .(
    chr, peak_start = start, peak_end = end
    #pref
  )]
  
  cat("  Generating pileups boxplot...\n")
  # make boxplot and select distance
  agg_both <- calculate.peak.pileup(nrwp, pu_both, logpileup=TRUE) 
  site_count <- get.site.count(nrwp,genome_sites)
  peak_sites <- merge(x = site_count, y = agg_both, by = c("chr", "peak_start", "peak_end"))
  sites_ctpile <- data.frame(count = merged$site.count, pileup = merged$pileup_score)
  box_plot_path <- paste(outpath,"/groupedbox_",len,".pdf",sep='')
  #plot.pileup.distributions(peak_sites$site.count, peak_sites$pileup.score,box_plot_path,chip_name=chip_name)
  
  cat("  Getting sites within distance for some peaks...\n")
  setDT(peak_sites, key = names(peak_sites))
  setkey(peak_sites, "chr", "peak_start","peak_end")

  for (ct in count_sites_per_peak){
    peaks_selected <- peak_sites[peak_sites$site_count == ct,] %>% select(-site_count) # count is not needed after this
    sites_all_dist <- get.consecutive.sites.in.peak(peaks_selected, genome_sites)
    sites_within_range <- get.probeseq.in.range(sites_all_dist, min_site_dist, max_site_dist, probe_len = probe_size, flank_size = probeseq_flank)
    sites_within_range$sites_in_peak = ct
    sites_within_range$peaklen = len
    idx <- idx + 1 # index starts with 1 in R
    datalist[[idx]] <- sites_within_range
  }
}
all_tbl <- rbindlist(datalist)  %>%
        mutate(s1 = s1 - seqstart + 1, 
               s2 = s2 - seqstart + 1,
               tf1 = tfname,
               tf2 = tfname)
all_tbl$key <- sprintf("sequence%s",seq.int(nrow(all_tbl)))
tbl_path <- paste(outpath,"/sites_all.tsv",sep='')
write.table(all_tbl,file=tbl_path,sep="\t",row.names = FALSE, quote = FALSE)
  