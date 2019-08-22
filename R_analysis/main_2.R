curwd <- "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"
setwd(curwd)

source("R_analysis/chip.info/R/cistreader.R")
source("R_analysis/chip.info/R/bsite.R")
source("R_analysis/chip.info/R/plotter.R")
source("R_analysis/chip.info/R/site_finder.R")

# ========= INPUT DIRECTORY ========= 
peaks_dir <- "../result/cistrome/peaks/hg19"
bed_path <- "../imads_preds/predictions/hg19_0005_Ets1_filtered.bed"
outpath <- "../result/cistrome"
tfname <- "cistrome_ets1"

# ========= INPUT PARAMETERS ========= 
probe_size <- 36
probeseq_flank <- 10 # n to the left and n to the right
spans <- c(50,100,150)
count_sites_per_peak <- c(2,3,4)
min_bsite_dist <- 1
max_bsite_dist <- 24

# ================== 

test.cistrome <- function(infile){
  cat(paste("Working on", infile, "\n"))
  tfid <- sapply(strsplit(infile,"_"), `[`, 1)
  tfname_id <- paste(tfname,tfid,sep="_")
  tf_inpath <- file.path(peaks_dir, infile)
  tf_outdir <- file.path(outpath, tfname_id, "analysis_result")
  sitefiles <- ""
  for (span in spans){
    cat(paste("  Iterating span", span, "...\n"))
    for (ct in count_sites_per_peak){
      cat(paste("    Iterating count", ct, "...\n"))
      sitefilename <- paste(tf_outdir,"/sites_within_d",ct,"_span",span,".tsv",sep='')
      sitefilepath <- paste(curwd,sitefilename,sep='/')
      sitefiles <- ifelse(sitefiles == "",sitefilepath,paste(sitefiles,sitefilepath,sep="\n"))
    }
  }
  return(sitefiles)
}

analyze.cistrome.peaks <- function(infile){
  cat(paste("Working on", infile, "\n"))
  tfid <- sapply(strsplit(infile,"_"), `[`, 1)
  tf_inpath <- file.path(peaks_dir, infile)
  tf_outdir <- file.path(outpath, tfid)
  
  # will create the directory if not exist
  dir.create(tf_outdir, showWarnings = FALSE)
  
  # first for every peak, we get the whole peak just to get the length distribution
  cat(paste("  Making the peak distances distribution\n"))
  peaklendist <- (read.cistrome.peak(tf_inpath) %>% mutate(len = end - start + 1))$len
  make.peaklen.dist.cistrome.plot(peaklendist, tf_outdir, chip_name=tfid)
  
  sitefiles <- ""
  for (span in spans){
    cat(paste("  Iterating span", span, "...\n"))
    # now we for every peak, we only get the span
    peak_df <- read.cistrome.peak(tf_inpath, span = span) 
    peak_bsites <- get.bsite.cist.count(peak_df, imads_df)
    box_plot_path <- paste(tf_outdir,"/groupedbox_",span,".pdf",sep='')
    make.intensity.dist.plot(peak_bsites, box_plot_path, tfid)
    
    setDT(peak_bsites, key = names(peak_bsites))
    setkey(peak_bsites, "chr", "peak.start","peak.end")
    for (ct in count_sites_per_peak){
      cat(paste("    Iterating count", ct, "...\n"))
      peaks_selected <- peak_bsites[peak_bsites$count == ct,] %>% select(-count) # count is not needed after this
      sites_all_dist <- get.bsites.in.peak(peaks_selected, imads_df, score_colname="intensity")
      sites_within_range <- get.bsites.within.range(sites_all_dist, min_bsite_distte_dist, max_bsite_dist, probeseq_flank)
      # add sequence id
      seqid <- paste("sequence",seq.int(nrow(sites_within_range)),sep="")
      sites_within_range <- cbind(seqid,sites_within_range)
      
      write.table(sites_all_dist,file=paste(tf_outdir,"/sites_all_d",ct,"_span",span,".tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
      sitefilename <- paste(tf_outdir,"/sites_within_d",ct,"_span",span,".tsv",sep='')
      write.table(sites_within_range,file=sitefilename,sep="\t",row.names = FALSE, quote = FALSE)
      
      sitefilepath <- paste(curwd,sitefilename,sep='/')
      sitefiles <- ifelse(sitefiles == "",sitefilepath,paste(sitefiles,sitefilepath,sep="\n"))
    }
  }
  return(sitefiles)
}

infiles <- list.files(path = peaks_dir, pattern = "\\.bed$")
#tfids <- sapply(strsplit(infiles,"_"), `[`, 1)

imads_df <- read.imads.bed(bed_path)

sitefiles <- sapply(infiles,test.cistrome)

write(sitefiles,file = paste(outpath,"/sitefiles_list.txt",sep=''))

