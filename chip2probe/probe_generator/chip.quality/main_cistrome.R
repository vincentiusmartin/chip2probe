# args <- commandArgs(trailingOnly=TRUE)
# print(args)
# setwd(args[1])
# cistrome_path <- args[2]
# imads_path <- args[3]
# outpath <- args[4]
# chip_name <- args[5]


setwd("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/probe_generator")
cistrome_path <- "../../resources/cistrome_files_ets1/37927.bed"
bed_path <- "../../resources/imads_files/predictions/hg19_0005_Ets1_filtered.bed"
outpath <- "result/cistrome_ets1_37927/analysis_result"
chip_name <- "cistrome_ets1_37927"


source("chip.quality/R/quality_check.R")

probe_size <- 36
probeseq_flank <- 10 # n to the left and n to the right
peaklen <- c(100,200,300)
count_sites_per_peak <- c(2,3,4)
min_bsite_dist <- 1
max_bsite_dist <- 24

cat("Reading binding sites prediction and input peaks...\n")
genome_sites <- read.sites.bed(bed_path)
nrwp_cistrome_all <- read.narrow.peak(cistrome_path, peak_type="all", show_peaklen=TRUE)
plot.peaklen.dist(list(nrwp_cistrome_all$peaklen), c("peaklen distribution"), paste(outpath,"/peaklen_dist.pdf",sep=''), chip_name)

datalist <- list()
idx <- 0
for (len in peaklen){
  cat("Working on analysis for peaklen",len,"...\n")
  nrwp <- read.narrow.peak(cistrome_path, peak_type="center", peak_length=len)
  peaklenlabel <- paste(chip_name,", peaklen ",len, sep='')
  site_count <- get.site.count(nrwp,genome_sites)
  cat("  Getting sites within distance for some peaks...\n")
  for (ct in count_sites_per_peak){
    peaks_selected <- site_count[site_count$site.count == ct,] %>% select(-site.count) # count is not needed after this
    sites_all_dist <- get.consecutive.sites.in.peak(peaks_selected, genome_sites)
    sites_within_range <- get.probeseq.in.range(sites_all_dist, min_bsite_dist, max_bsite_dist, probe_len = probe_size, flank_size = probeseq_flank)
    sites_within_range$sites.in.peak <- ct
    sites_within_range$peaklen <- len
    idx <- idx + 1 # index starts with 1 in R
    datalist[[idx]] <- sites_within_range
  }
}
all_tbl <- rbindlist(datalist)
all_tbl$key <- sprintf("sequence%s",seq.int(nrow(all_tbl)))
tbl_path <- paste(outpath,"/sites_all.tsv",sep='')
write.table(all_tbl,file=tbl_path,sep="\t",row.names = FALSE, quote = FALSE)
