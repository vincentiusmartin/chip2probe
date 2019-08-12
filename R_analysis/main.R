#.libPaths( c( .libPaths(), "/data/gordanlab/vincentius/cooperative_probe/packages/Rlib") )
#curwd <- "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"
#setwd(curwd)
source("R_analysis/chip.info/R/chipreader.R")
source("R_analysis/chip.info/R/bsite.R")
source("R_analysis/chip.info/R/plotter.R")
source("R_analysis/chip.info/R/site_finder.R")

#args = commandArgs(trailingOnly=TRUE)
#print(args)
#setwd("/Users/vincentiusmartin/Research/chip2gcPBM/src/../")


# #setwd("/Users/vincentiusmartin/Research/chip2gcPBM/src")
# #curwd <- "/Users/vincentiusmartin/Research/chip2gcPBM/src"
# pu1_path <- "../result/ets1_k562/macs_result/ets1_k562_r1_treat_pileup.bdg"
# pu2_path <- "../result/ets1_k562/macs_result/ets1_k562_r2_treat_pileup.bdg"
# pu_both_path <- "../result/ets1_k562/macs_result/ets1_k562_bothrs_treat_pileup.bdg"
# nrwp_preidr_path <- "../result/ets1_k562/macs_result/ets1_k562_bothrs_peaks.narrowPeak"
# nrwp_postidr_path <- "../result/ets1_k562/idr_result/idr_001p_wlist.005i"
# bed_path <- "../imads_files/predictions/hg19_0005_Ets1_filtered.bed"
# outpath <- "../result/ets1_A549/analysis_result"
# chip_name <- "ets1_A549"

args = commandArgs(trailingOnly=TRUE)
print(args)
setwd(args[1])
curwd <- args[1]
pu1_path <- args[2]
pu2_path <- args[3]
pu_both_path <- args[4]
nrwp_preidr_path <- args[5]
nrwp_postidr_path <- args[6]
bed_path <- args[7]
outpath <- args[8]
chip_name <- args[9]

probe_size <- 36
probeseq_flank <- 10 # n to the left and n to the right
spans <- c(50,100,150)
count_sites_per_peak <- c(2,3,4)
min_bsite_dist <- 1
max_bsite_dist <- 24

# ===== END OF CONFIGURATION ===== 

cat("Reading binding sites prediction...\n")
imads_bsite <- read.imads.bed(bed_path)
cat("Reading pileup files...\n")
pu1 <- read.pileup(pu1_path)
pu2 <- read.pileup(pu2_path)
pu_both <- read.pileup(pu_both_path)

nrwp_preidr <- read.narrow.peak(nrwp_preidr_path)
nrwp_postidr <- read.narrow.peak(nrwp_postidr_path)
wdist_preidr <- nrwp_preidr %>% mutate(len = end - start + 1)
wdist_postidr <- nrwp_postidr %>% mutate(len = end - start + 1)
write.table(wdist_preidr,file=paste(outpath,"/narrowpeak_preidr.tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
write.table(wdist_postidr,file=paste(outpath,"/narrowpeak_postidr.tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
write.peak.info(wdist_preidr$len, wdist_postidr$len, paste(outpath,"/peakinfo.txt",sep=''))

make.peaklen.dist.plot(wdist_preidr$len, wdist_postidr$len, chip_name, outpath)

sitefiles <- ""
for (span in spans){
  cat("Working on analysis for span",span,"...\n")
  nrwp <- read.narrow.peak(nrwp_postidr_path, span=span)
  spanlabel <- paste(chip_name,", span ",span, sep='')
  
  cat("  Generating correlation plot...\n")
  agg1 <- get.pileups(nrwp, pu1, logpileup=FALSE)
  agg2 <- get.pileups(nrwp, pu2, logpileup=FALSE)
  merged <- merge(x = agg1, y = agg2, by = c("chr", "peak.start", "peak.end"))
  write.table(merged,file=paste(outpath,"/pileup_scores_span",span,".tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
  
  corr_plot_path <- paste(outpath,"/corr_plot_",span,".pdf",sep='')
  # make the correlation plot using the logged value
  make.corr.plot(log(merged$pileups.x+1), log(merged$pileups.y+1), corr_plot_path, spanlabel)
  
  cat("  Generating pileups boxplot...\n")
  # make boxplot and select distance
  agg_both <- get.pileups(nrwp, pu_both, logpileup=TRUE) 
  bsite_count <- get.bsite.count(nrwp,imads_bsite)
  peak_bsites <- merge(x = bsite_count, y = agg_both, by = c("chr", "peak.start", "peak.end"))
  bsites_ctpile <- peak_bsites[ , c("count","pileups")]
  
  box_plot_path <- paste(outpath,"/groupedbox_",span,".pdf",sep='')
  make.pileup.dist.plot(bsites_ctpile, box_plot_path, spanlabel)

  cat("  Getting sites within distance for some peaks...\n")
  setDT(peak_bsites, key = names(peak_bsites))
  setkey(peak_bsites, "chr", "peak.start","peak.end")
  for (ct in count_sites_per_peak){
    # Get binding sites in a peak 
    peaks_selected <- peak_bsites[peak_bsites$count == ct,] %>% select(-count) # count is not needed after this
    sites_all_dist <- get.bsites.in.peak(peaks_selected, imads_bsite)
    sites_within_range <- get.bsites.within.range(sites_all_dist, min_bsite_dist, max_bsite_dist, probeseq_flank)

    # using http to get genomic seq: https://bioinformatics.stackexchange.com/questions/2543/way-to-get-genomic-sequences-at-given-coordinates-without-downloading-fasta-file
    write.table(sites_all_dist,file=paste(outpath,"/sites_all_d",ct,"_span",span,".tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
    #sitefilename <- paste(outpath,"/sites_within_d",ct,"_span",span,".tsv",sep='')
    write.table(sites_within_range,file=sitefilename,sep="\t",row.names = FALSE, quote = FALSE)
    
    #sitefilepath <- paste(curwd,sitefilename,sep='/')
    #sitefiles <- ifelse(sitefiles == "",sitefilename,paste(sitefilepath,sitefilename,sep="\n"))
  }
}

sitefiles <- ""
for (span in spans){
  for (ct in count_sites_per_peak){
    sitefilename <- paste(outpath,"/sites_within_d",ct,"_span",span,".tsv",sep='')
    sitefilepath <- paste(curwd,sitefilename,sep='/')
    print(sitefilepath)
    sitefiles <- ifelse(sitefiles == "",sitefilepath,paste(sitefiles,sitefilepath,sep="\n"))
  }
}
# write all file with sequences to a file to enable parsing
write(sitefiles,file = paste(outpath,"/sitefiles_list.txt",sep=''))
