library(data.table)
library(dplyr)
library(ggplot2)
library("BSgenome.Hsapiens.UCSC.hg19")

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
#' read.imads.bed()
calculate.peak.pileup <- function(nrwp_df, pu_df, logpileup=FALSE, jointype="any"){
  pu_nrwp <- foverlaps(pu_df, nrwp_df, nomatch=NULL, type=jointype)[, .(
    chr, peak.start = start, peak.end = end,
    pileup.start = i.start, pileup.end = i.end,
    pileup
  )] %>% 
    mutate(
      # Handle edge cases here, consider partial reads by bordering with both ends of the peak
      pileup.start = ifelse(pileup.start < peak.start, peak.start, pileup.start),
      pileup.end = ifelse(pileup.end > peak.end, peak.end, pileup.end),
      pileup.score = (pileup.end - pileup.start) * pileup
    ) %>%
    select(-c(pileup.end,pileup.start,pileup))
  
  setDT(pu_nrwp, key = c("chr", "peak.start","peak.end"))
  agg <- pu_nrwp[, lapply(.SD, sum), by=list(chr,peak.start,peak.end)]
  if(logpileup){
    # +1 to avoid -inf on log
    agg$pileup.score = log(agg$pileup.score+1)
  }
  
  return(agg)
}

# ================= iMADS part ================= 

#' Read imads bed
#'
#' This function reads a narrow peak file and make a data table from it.
#' @param bed_path The input imads bed path prediction.
#' @keywords 
#' @export
#' @examples
#' read.imads.bed()
read.imads.bed <- function(bed_path){
  imads_sites <- fread(bed_path,sep="\t",header=FALSE)[,1:4]
  # the fifth column is basically the same, so just use 4 columns
  colnames(imads_sites) <- c("chr","start","end","pref")
  setDT(imads_sites, key = c("chr", "start","end"))
  return(imads_sites)
}

# ================= Analysis with imads and macs file =================


#' Get binding site count
#'
#' For each peak, gets the number of binding site--predicted by imads.
#' @param peaklen_lists The input imads bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.site.count()
get.site.count <- function(nrwp_df, imads_df){
  sites_peak <- foverlaps(imads_df,nrwp_df, nomatch=NULL)[, .(
    chr, peak.start = start, peak.end = end,
    pref
  )]
  
  sites_peak$count <- ifelse(is.na(sites_peak$pref), 0, 1)
  sites_peak <- sites_peak %>% select(-pref)
  counted <- sites_peak[, lapply(.SD, sum), by=list(chr,peak.start,peak.end)]
  colnames(counted) <- c("chr","peak.start","peak.end","site.count")
  
  return(counted)
}

#' consecutive sites in peak
#'
#' For each peak, gets all binding sites available
#' @param pileups_df The input imads bed path prediction.
#' @param imads_df The input imads bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.bsites.in.peak()
get.sites.in.peak <- function(peak_df, imads_df, join_type="within"){
  peaks <- peak_df
  sites <- imads_df
  colnames(peaks)[c(2,3)] <- c("peak.start","peak.end")
  colnames(sites)[c(2,3)] <- c("bsite.start","bsite.end")
  sites_in_peak <- foverlaps(sites, peaks, nomatch=NULL, type=join_type)
  return(sites_in_peak)
}

#' consecutive sites in peak
#'
#' For each peak, gets all binding sites available
#' @param pileups_df The input imads bed path prediction.
#' @param imads_df The input imads bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.bsites.in.peak()
get.consecutive.sites.in.peak <- function(peak_df, imads_df, score_colname = "pileups"){
  peaks <- peak_df
  sites <- imads_df
  names(peaks)[c(2,3)] <- c("peak.start","peak.end")
  names(sites)[c(2,3)] <- c("bsite.start","bsite.end")
  setDT(peaks, key=c("chr","peak.start","peak.end"))
  setDT(sites, key=c("chr","bsite.start","bsite.end"))
  sites_selected <- foverlaps(sites, peaks, nomatch=NULL, type="within")
  # order the chromosome
  chr_order<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  sites_selected$chr <- factor(sites_selected$chr, levels=chr_order)
  sites_selected <- sites_selected[order(sites_selected$chr),]
  
  sites_in_peak_df <- sites_selected %>%
    group_by(chr, peak.start, peak.end) %>%
    mutate(distance = (bsite.start + bsite.end) / 2 - lag((bsite.start + bsite.end)/2),
           bsite.start.1 = lag(bsite.start),
           bsite.end.1 = lag(bsite.end),
           bsite.start.2 = bsite.start,
           bsite.end.2 = bsite.end,
           bsite.pref.1 = lag(pref),
           bsite.pref.2 = pref) %>%
    filter(!is.na(distance)) %>%
    select(-c(bsite.start, bsite.end, pref)) %>%
    select(-distance,everything()) #%>% # just to put distance as the last col
  
  return(sites_in_peak_df)
}


#' Get binding sites in peak
#'
#' Important columns needed are:
#' 
#' 
#' @param sites_in_peak_df The input imads bed path prediction.
#' @param min_dist The input imads bed path prediction.
#' @param max_dist The input imads bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.bsites.in.peak()
get.probeseq.in.range <- function(sites_in_peak_df, min_dist, max_dist, probe_len = 36, flank_size = 0){
  # TODO: filter chrM?
  # Get binding sites in min_dist .. max_dist
  # if no sites are within distance, a warning message would show up:
  # Factor `chr` contains implicit NA, consider using `forcats::fct_explicit_na`
  sites_within_range <- sites_in_peak_df %>%
    filter(distance >= min_dist & distance <= max_dist)  %>%
    # rowwise() %>% # need to calculate each row separately <- DON'T USE THIS
    mutate(
      s1 = ifelse(bsite.start.1 < bsite.start.2, ((bsite.start.1 + bsite.end.1) / 2), ((bsite.start.2 + bsite.end.2) / 2)),
      s2 = ifelse(bsite.start.1 > bsite.start.2, ((bsite.start.1 + bsite.end.1) / 2), ((bsite.start.2 + bsite.end.2) / 2)),
      seqstart = s1 - ceiling((probe_len - distance) / 2) + 1,
      seqend =  s2 + floor((probe_len - distance) / 2)
    ) %>%
    select(-c(s1,s2)) # made just as a tmp variable so remove here
  
  # don't add +1 on everything to conform with one index
  sequences <- getSeq(
    Hsapiens,
    sites_within_range$chr,
    start = sites_within_range$seqstart,
    end = sites_within_range$seqend
  )
  flank_left <- getSeq(
    Hsapiens,
    sites_within_range$chr,
    start = sites_within_range$seqstart - flank_size,
    end = sites_within_range$seqstart - 1
  )
  flank_right <- getSeq(
    Hsapiens,
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
#' @param sites_in_peak_df The input imads bed path prediction.
#' @param min_dist The input imads bed path prediction.
#' @param max_dist The input imads bed path prediction.
#' @keywords 
#' @export
#' @examples
#' get.bsites.in.peak()
get.closest.site.in <- function(site_df1, site_df2){
  if (!all(c("chr", "bsite.start", "bsite.end") %in% colnames(site_df1)) && 
      !all(c("chr", "start","end") %in% colnames(site_df1))
  ){
    stop("Could not find columns with genomic coordinates in the first table.") 
  }
  
  if (!all(c("chr", "bsite.start", "bsite.end") %in% colnames(site_df2)) && 
      !all(c("chr", "start","end") %in% colnames(site_df2))
  ){
    stop("Could not find columns with genomic coordinates in the second table.") 
  }
  
  df1 <- copy(site_df1)
  df2 <- copy(site_df2)
  
  if (all(c("chr", "start","end") %in% colnames(df1)) && 
      !all(c("chr", "bsite.start", "bsite.end") %in% colnames(df1))){
    setnames(df1, old = c("start", "end"), new = c("bsite.start", "bsite.end"))
  }
  if (all(c("chr", "start","end") %in% colnames(df2)) && 
      !all(c("chr", "bsite.start", "bsite.end") %in% colnames(df2))){
    setnames(df2, old = c("start", "end"), new = c("bsite.start", "bsite.end"))
  }
  
  same_cname <- intersect(colnames(df1),colnames(df2))
  same_cname <- same_cname[!same_cname %in% "chr"]
  cname1 <- paste(same_cname, "1", sep=".")
  cname2 <- paste(same_cname, "2", sep=".")
  setnames(df1, old = same_cname, new = cname1)
  setnames(df2, old = same_cname, new = cname2)
  
  df1 <- df1 %>% mutate(
    pos = (bsite.start.1 + bsite.end.1) %/% 2,
    pos.1 = pos
  )
  df2 <- df2 %>% mutate(
    pos = (bsite.start.2 + bsite.end.2) %/% 2,
    pos.2 = pos
  )
  
  setDT(df1, key = c("chr", "bsite.start.1","bsite.end.1"))
  setDT(df2, key = c("chr", "bsite.start.2","bsite.end.2"))
  
  # rolling merge will be done by matching chr and then do rolling join 
  # using st (i.e. "site") that's why the join is on chr,st
  rolled <- df2[df1, roll = "nearest", on = .(chr,pos)] %>% 
    mutate(distance = abs(pos.1-pos.2)) %>%
    select(-c("pos"))
  setDT(rolled)
  
  return(rolled)
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
  bsites_ctpile <- data.frame(count = count_vec, pileup = pileup_vec)
  merged_copy <- bsites_ctpile
  merged_copy$count <- "all"
  merged_duplicated <- rbind(bsites_ctpile,merged_copy)
  countvec <- c("all",as.character(sort(unique(bsites_ctpile$count)))) # need this to order xlabel
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
