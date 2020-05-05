library(data.table)
library(ggplot2)
library(dplyr)
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")
source("/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/chip2probe/probe_generator/chip.quality/R/kompas.R")

# ================= ChIP-seq part ================= 

#' Read a narrow peak file
#'
#' This function reads a narrow peak file and make a data table from it.
#' @param nrwp_path The input narrow peak file path.
#' @param peak_type Default value is \code{all}. Allowed values are \code{all},
#'        \code{summit}, or \code{center}.
#'        With \code{all} it will use the provided peak start and peak end from 
#'        the input file. If \code{center} is used then the peak start and end
#'        will be around the midpoint of the peak. For \code{summit}, the 
#'        function makes peak_start and end around the original peak start + 
#'        \code{summit_col_idx}.
#' @param peak_length The length of the peak when \code{peak_type} is not set to 
#'        all.
#' @param summit_col_idx When \code{peak_type == summit}, this denotes the column
#'        with the summit position. If using macs2, this will be the 10th column,
#'        and thus \code{10} is the default value.
#' @param show_peaklen If \code{TRUE} add a peaklen column. Useful when 
#'        \code{peak_type == all} to get the length of each peak.
#' @keywords narrow peak, data table
#' @export
#' @examples
#' read.narrow.peak()
read.narrow.peak <- function(nrwp_path, peak_type = "all", peak_length = -1, summit_col_idx = 10, 
                             show_peaklen = FALSE){
  # Error checking for the parameter
  if(!peak_type %in% c("all", "summit", "center")){
    stop("peak_type needs to be either \"all\", \"summit\", or \"center\"") 
  }
  if(peak_type != "all" && peak_length == -1){
    stop("peak length needs to be set if the peak_type is not \"all\"") 
  }
  if(peak_type == "all" && peak_length != -1){
    warnings(cat("Peak length is set to",peak_length,",but the peak type is all. Please change the peak type to summit or center.")) 
  }
  
  span <- peak_length %/% 2
  if(peak_type == "summit"){
    nrwp <- fread(nrwp_path, select = c(1,2,3,5,summit_col_idx))
    colnames(nrwp) <- c("chr","start","end","score","summit")
    summit_pos <- nrwp$start + nrwp$summit
    nrwp$start <- summit_pos - span
    nrwp$end <- summit_pos + span - 1
    # drop the summit column since it is no longer needed
    nrwp <- subset(nrwp, select = -c(summit)) 
  }else{ # all or center
    nrwp <- fread(nrwp_path, select = c(1,2,3,5))
    colnames(nrwp) <- c("chr","start","end","score")
    if(peak_type == "center"){
      # some peak files like from cistrome don't have peak summit,
      # so we can use the peak center.
      peak_center <- (nrwp$start + nrwp$end) %/% 2
      nrwp$start <- peak_center - span
      nrwp$end <- peak_center + span - 1
    }
  }
  
  allowed_chr = c(paste("chr",1:22,sep=""),"chrX","chrY")
  nrwp <- filter(nrwp, chr %in% allowed_chr) #filter chrM and other unsupported chr
  #nrwp <- nrwp[(nrwp$chr != "chrM"),]
  if(show_peaklen){
    nrwp <- nrwp %>% mutate(peaklen = end - start + 1) # +1 because the peaklen should be inclusive
  }
  setDT(nrwp, key = c("chr", "start","end"))
  return(nrwp)
}

#' Read pileup file
#'
#' This function reads a pileup file
#' @param pileup_path The input pileup file.
#' @keywords 
#' @export
#' @examples
#' read.pileup()
read.pileup <- function(pileup_path){
  pu <- fread(pileup_path)
  colnames(pu) <- c("chr","start","end","pileup")
  setDT(pu, key = c("chr", "start","end"))
  return(pu)
}

#' Calculate peak pileup
#'
#' Calculates the pileup value from every peak
#' @param nrwp_df explain
#' @param pu_df explain.
#' @param logpileup explain.
#' @param jointype explain.
#' @keywords 
#' @export
#' @examples
#' calculate.peak.pileup()
calculate.peak.pileup <- function(nrwp_df, pu_df, logpileup=FALSE, jointype="any"){
  pu_nrwp <- foverlaps(pu_df, nrwp_df, nomatch=NULL, type=jointype)[, .(
    chr, peak_start = start, peak_end = end,
    pileup_start = i.start, pileup_end = i.end,
    pileup
  )] %>% 
    mutate(
      # Handle edge cases here, consider partial reads by bordering with both ends of the peak
      pileup_start = ifelse(pileup_start < peak_start, peak_start, pileup_start),
      pileup_end = ifelse(pileup_end > peak_end, peak_end, pileup_end),
      pileup_score = (pileup_end - pileup_start) * pileup
    ) %>%
    select(-c(pileup_end,pileup_start,pileup))
  
  setDT(pu_nrwp, key = c("chr", "peak_start","peak_end"))
  agg <- pu_nrwp[, lapply(.SD, sum), by=list(chr,peak_start,peak_end)]
  if(logpileup){
    # +1 to avoid -inf on log
    agg$pileup_score = log(agg$pileup_score+1)
  }
  
  return(agg)
}

# ================= sites bed file part ================= 

#' Read sites bed
#'
#' This function reads a narrow peak file and make a data table from it.
#' @param bed_path The input sites bed path prediction.
#' @keywords 
#' @export
#' @examples
#' read.imads.bed()
read.sites.bed <- function(bed_path){
  genome_sites <- fread(bed_path,sep="\t",header=FALSE)[,1:3] # imads: [,1:4]
  # the fifth column is basically the same, so just use 4 columns
  colnames(genome_sites) <- c("chr","start","end")
  setDT(genome_sites, key = c("chr", "start","end"))
  return(genome_sites)
}

# ================= Read sites bed using configuration ================= 
#' peak file is only needed when kompas is used
get.genome.sites <- function(config, tfname, peak_file="", genomever="hg19"){
  sites_type <- config$sitecall_mode
  
  if(sites_type == "imads"){
    cat("Using iMADS bed file", tfname ,"...\n")
    bedpath <- config$imads[[tfname]]$bedpath
    genome_sites <- read.sites.bed(bedpath)
  }else{
    align_path <- config$kompas[[tfname]]$alignpath
    pwm_start <- config$kompas[[tfname]]$pwm_core_start
    pwm_end <- config$kompas[[tfname]]$pwm_core_end
    pwm_pos <- c(pwm_start, pwm_end)
    cat("Running Kompas for ", tfname,"...\n",sep="")
    genome_sites <- run.kompas(align_path, peak_file, pwm_pos ,0.4, genomever=genomever)
    setDT(genome_sites, key = c("chr", "start", "end"))
  }
  
  return(genome_sites)
}

# ================= Analysis with genome sites and macs file =================

#' Get binding site count
#'
#' For each peak, gets the number of binding site--predicted by imads.
#' @param peaklen_lists The input imads bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.site.count()
get.site.count <- function(nrwp_df, sites_df){ #here todo
  sites_peak <- foverlaps(sites_df,nrwp_df, nomatch=NULL)[, .(
    chr, peak_start = start, peak_end = end
    #pref
  )]
  
  sites_peak$count <- 1 #ifelse(is.na(sites_peak$pref), 0, 1)
  #sites_peak <- sites_peak %>% select(-pref)
  counted <- sites_peak[, lapply(.SD, sum), by=list(chr,peak_start,peak_end)]
  colnames(counted) <- c("chr","peak_start","peak_end","site_count")
  
  return(counted)
}

#' consecutive sites in peak
#'
#' For each peak, gets all binding sites available
#' @param pileups_df The input sites bed path prediction.
#' @param imads_df The input sites bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.sites.in.peak()
get.sites.in.peak <- function(peak_df, imads_df, join_type="within"){
  peaks <- peak_df
  sites <- imads_df
  colnames(peaks)[c(2,3)] <- c("peak_start","peak_end")
  colnames(sites)[c(2,3)] <- c("site_start","site_end")
  sites_in_peak <- foverlaps(sites, peaks, nomatch=NULL, type=join_type)
  return(sites_in_peak)
}

#' consecutive sites in peak
#'
#' For each peak, gets all binding sites available
#' @param pileups_df The input sites bed path prediction.
#' @param imads_df The input sites bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.sites.in.peak()
get.consecutive.sites.in.peak <- function(peak_df, sites_df, score_colname = "pileups"){
  peaks <- peak_df
  sites <- sites_df
  names(peaks)[c(2,3)] <- c("peak_start","peak_end")
  names(sites)[c(2,3)] <- c("site_start","site_end")
  setDT(peaks, key=c("chr","peak_start","peak_end"))
  setDT(sites, key=c("chr","site_start","site_end"))
  sites_selected <- foverlaps(sites, peaks, nomatch=NULL, type="within")
  # order the chromosome
  chr_order<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  sites_selected$chr <- factor(sites_selected$chr, levels=chr_order)
  sites_selected <- sites_selected[order(sites_selected$chr),]
  
  sites_in_peak_df <- sites_selected %>%
    group_by(chr, peak_start, peak_end) %>%
    mutate(distance = (site_start + site_end) / 2 - lag((site_start + site_end)/2),
           site_start_1 = lag(site_start),
           site_end_1 = lag(site_end),
           site_start_2 = site_start,
           site_end_2 = site_end
           #site.pref.1 = lag(pref),
           #site.pref.2 = pref
           ) %>%
    filter(!is.na(distance)) %>%
    select(-c(site_start, site_end)) %>%
    select(-distance,everything()) #%>% # just to put distance as the last col
  
  return(sites_in_peak_df)
}


#' Get binding sites in peak
#'
#' Important columns needed are:
#' 
#' 
#' @param sites_in_peak_df The input sites bed path prediction.
#' @param min_dist The input sites bed path prediction.
#' @param max_dist The input sites bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.probeseq.in.range 
get.probeseq.in.range <- function(sites_in_peak_df, min_dist, max_dist, 
                                  probe_len = 36, flank_size = 0, genomever="hg19"){
  # TODO: filter chrM?
  # Get binding sites in min_dist .. max_dist
  # if no sites are within distance, a warning message would show up:
  # Factor `chr` contains implicit NA, consider using `forcats::fct_explicit_na`
  sites_within_range <- sites_in_peak_df %>%
    filter(distance >= min_dist & distance <= max_dist)  %>%
    # rowwise() %>% # need to calculate each row separately <- DON'T USE THIS
    mutate(
      mid_1 = ifelse(site_start_1 < site_start_2, ((site_start_1 + site_end_1) %/% 2), ((site_start_2 + site_end_2) %/% 2)),
      mid_2 = ifelse(site_start_1 > site_start_2, ((site_start_1 + site_end_1) %/% 2), ((site_start_2 + site_end_2) %/% 2)),
      seqstart = mid_1 - ceiling((probe_len - distance) / 2) + 1,
      seqend =  mid_2 + floor((probe_len - distance) / 2)
    )
  
  if (genomever == "hg19") {
    hs_gene <- BSgenome.Hsapiens.UCSC.hg19
  }else{
    hs_gene <- BSgenome.Hsapiens.UCSC.hg38
  }
  # don't add +1 on everything to conform with one index
  sequences <- getSeq(
    hs_gene,
    sites_within_range$chr,
    start = sites_within_range$seqstart,
    end = sites_within_range$seqend
  )
  
  flank_left <- getSeq(
    hs_gene,
    sites_within_range$chr,
    start = sites_within_range$seqstart - flank_size,
    end = sites_within_range$seqstart - 1
  )
  flank_right <- getSeq(
    hs_gene,
    sites_within_range$chr,
    start = sites_within_range$seqend,
    end = sites_within_range$seqend + flank_size - 1
  )
  
  sites_within_range$sequence <- as.character(sequences)
  sites_within_range$flank_left <- as.character(flank_left)
  sites_within_range$flank_right <- as.character(flank_right)
  
  return(sites_within_range)
}

# ================= Heterotypic cluster =================

#' Get binding sites in peak
#'
#' Important columns needed are:
#' 
#' 
#' @param sites_in_peak_df The input sites bed path prediction.
#' @param min_dist The input sites bed path prediction.
#' @param max_dist The input sites bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.sites.in.peak()
get.closest.site.in <- function(site_df1, site_df2){
  if (!all(c("chr", "site.start", "site.end") %in% colnames(site_df1)) && 
      !all(c("chr", "start","end") %in% colnames(site_df1))
  ){
    stop("Could not find columns with genomic coordinates in the first table.") 
  }
  
  if (!all(c("chr", "site.start", "site.end") %in% colnames(site_df2)) && 
      !all(c("chr", "start","end") %in% colnames(site_df2))
  ){
    stop("Could not find columns with genomic coordinates in the second table.") 
  }
  
  df1 <- copy(site_df1)
  df2 <- copy(site_df2)
  
  if (all(c("chr", "start","end") %in% colnames(df1)) && 
      !all(c("chr", "site_start", "site_end") %in% colnames(df1))){
    setnames(df1, old = c("start", "end"), new = c("site_start", "site_end"))
  }
  if (all(c("chr", "start","end") %in% colnames(df2)) && 
      !all(c("chr", "site_start", "site_end") %in% colnames(df2))){
    setnames(df2, old = c("start", "end"), new = c("site_start", "site_end"))
  }
  
  # get which columns to use for join
  same_cname <- intersect(colnames(df1),colnames(df2))
  same_cname <- same_cname[!same_cname %in% "chr"]
  cname1 <- paste(same_cname, "1", sep=".")
  cname2 <- paste(same_cname, "2", sep=".")
  setnames(df1, old = same_cname, new = cname1)
  setnames(df2, old = same_cname, new = cname2)
  
  # for this, assume that site position is in the middle
  df1 <- df1 %>% mutate(
    pos = (site_start_1 + site_end_1) %/% 2,
    pos_1 = pos
  )
  df2 <- df2 %>% mutate(
    pos = (site_start_2 + site_end_2) %/% 2,
    pos_2 = pos
  )
  
  setDT(df1, key = c("chr", "site_start_1","site_end_1"))
  setDT(df2, key = c("chr", "site.start_2","site_end_2"))
  
  # rolling merge will be done by matching chr and then do rolling join 
  # using st (i.e. "site") that's why the join is on chr,st
  rolled <- df2[df1, roll = "nearest", on = .(chr,pos)] %>% 
    mutate(distance = abs(pos.1-pos.2)) %>%
    select(-c("pos"))
  setDT(rolled)
  
  return(rolled)
}

#' Get sites with overlapping peaks from two tf tables
#'
#' @param peaklen_lists the list of peak length. The input MUST BE a list where
#'        each element is a vector of peak length.
#' @param upper_coord_limit can set to -1 if desired.
#' @keywords 
#' @export
#' @examples
#' get.sites.with.overlapping.peaks()
get.heterotypic.sites <- function(peaks1, peaks2, sites1, sites2, dist_range){
  peak_all_intersect <- foverlaps(peaks1, peaks2, nomatch=0)[, .(
    chr, peak_start_1 = i.peak_start, peak_end_1 = i.peak_end,
    peak_start_2 = peak_start, peak_end_2 = peak_end
  )]
  # get the list of peaks that overlap in both
  peak1_intersect <- peak_all_intersect[, c("chr", "peak_start_1", "peak_end_1")]
  setDT(peak1_intersect, key=c("chr", "peak_start_1", "peak_end_1"))
  peak2_intersect <- peak_all_intersect[, c("chr", "peak_start_2", "peak_end_2")]
  setDT(peak2_intersect, key=c("chr", "peak_start_2", "peak_end_2"))
  
  peak1_sites <- foverlaps(sites1, peak1_intersect, nomatch=0, type="within")[, .(
    chr, peak_start_1, peak_end_1,
    site_start_1 = start, site_end_1 = end
  )]
  peak1_sites <- peak1_sites %>% mutate(
    # need to make mid column for rolling join and mid.1 for distance calculation
    mid = (site_start_1 + site_end_1) %/% 2,
    mid_1 = mid
  )
  setDT(peak1_sites, key=c("chr", "peak_start_1", "peak_end_1"))
  
  peak2_sites <- foverlaps(sites2, peak2_intersect, nomatch=0, type="within")[, .(
    chr, peak_start_2, peak_end_2,
    site_start_2 = start, site_end_2 = end
  )]
  peak2_sites <- peak2_sites %>% mutate(
    mid = (site_start_2 + site_end_2) %/% 2,
    mid_2 = mid
  )
  setDT(peak2_sites, key=c("chr", "peak_start_2", "peak_end_2"))
  
  # get sites close to tf1
  sites_close_to_1 <- peak1_sites[peak2_sites, roll = "nearest", on = .(chr,mid)] %>% 
    select(-c("mid")) %>%
    mutate(distance = abs(mid_1-mid_2)) %>%
    filter(distance > dist_range[1] & distance < dist_range[2])
  sites_close_to_2 <- peak2_sites[peak1_sites, roll = "nearest", on = .(chr,mid)] %>% 
    select(-c("mid")) %>%
    mutate(distance = abs(mid_1-mid_2)) %>%
    filter(distance > dist_range[1] & distance < dist_range[2])
  
  # merge both tables, take unique values
  all_sites <- rbind(sites_close_to_1, sites_close_to_2) %>% unique()
  
  chr_order<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  all_sites$chr <- factor(all_sites$chr, levels=chr_order)
  all_sites <- all_sites  %>% arrange(chr, peak_start_1, peak_end_1)
}

# ================= Plotting part =================

#' Make a plot of peak length distribution
#'
#' @param peaklen_lists the list of peak length. The input MUST BE a list where
#'        each element is a vector of peak length.
#' @param upper_coord_limit can set to -1 if desired.
#' @keywords 
#' @export
#' @examples
#' make.peaklen.dist.plot()
plot.peaklen.dist <- function(peaklen_lists, vec_labels, outpath, chipname = "", upper_coord_limit = 1500){
  upper_coord <- ifelse(upper_coord_limit != -1, min(c(sapply(peaklen_lists, max),upper_coord_limit)), sapply(peaklen_lists, max))
  pl_dist <- stack(setNames(peaklen_lists, vec_labels))
  names(pl_dist) <- c("values", "legend")
  peaklen_dist_plot <- ggplot(data = pl_dist, aes(x=values,colour=legend,fill=legend)) +
    #geom_histogram(binwidth=50, size=0.1, alpha=0.4) +
    geom_density(alpha=0.45)  + # overlay with transparent density plot
    coord_cartesian(xlim=c(min(unlist(peaklen_lists)),upper_coord)) +
    labs(title="ChIP-seq peak length distribution",
         subtitle=chipname,
         x="peak length") +
    theme(legend.position=c(.8,.8), legend.title=element_blank(), legend.spacing.x=unit(0.2,'cm'))
  ggsave(outpath)
}

#' Make correlation plot
#'
#' Make correlation plot between 2 ChIP-seq replicates
#' @param xvals explain
#' @param yvals explain.
#' @param path explain.
#' @param chipname explain.
#' @keywords 
#' @export
#' @examples
#' make.corr.plot()
plot.corr <- function(xvals, yvals, path, chip_name = ""){
  pu_cor <- cor(xvals,yvals)
  corr_plot <- ggplot() +
    geom_point(aes(x=xvals, y=yvals), size=0.5, colour="blue") +
    geom_abline(colour="red") +
    annotate("text", label = sprintf("R² = %.2f",pu_cor^2), x = min(xvals) + 0.5, y = max(yvals))  +
    labs(title="Correlation between replicates",
         subtitle=chip_name,
         y="Replicate 2",
         x="Replicate 1")
  ggsave(path)
}

#https://stackoverflow.com/questions/28846348/add-number-of-observations-per-group-in-ggplot2-boxplot
count.n <- function(x) {
  return(c(y = median(x)*1.03, label = length(x)))
}

#' Make pileup distribution plot
#'
#' Make pileup distribution plot with count for each pileup
#' @param xvals explain
#' @param yvals explain.
#' @param path explain.
#' @param chipname explain.
#' @keywords 
#' @export
#' @examples
#' make.pileup.dist.plot()
plot.pileup.distributions <- function(count_vec, pileup_vec, outpath, chip_name="", ylabel="log(pileup)"){
  # make boxplot, since we need pileup for the whole distribution, make a copy of the whole table
  # as group "all"
  sites_ctpile <- data.frame(count = count_vec, pileup = pileup_vec)
  merged_copy <- sites_ctpile
  merged_copy$count <- "all"
  merged_duplicated <- rbind(sites_ctpile,merged_copy)
  countvec <- c("all",as.character(sort(unique(sites_ctpile$count)))) # need this to order xlabel
  merged_duplicated$count <- factor(merged_duplicated$count, levels = countvec)
  pudist_plot <- ggplot(merged_duplicated, aes(x=count, y=pileup, fill=count)) +
    geom_boxplot() +
    stat_summary(
      # count.n is a function defined above to adjust the position of the count label
      fun.data = count.n,
      geom = "text"
    ) +
    labs(title="Pileup distributions",
         subtitle=chip_name,
         y=ylabel,
         x="#sites in a peak") +
    theme(legend.position = "none")
  ggsave(outpath)
}