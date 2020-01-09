library(data.table)
library(dplyr)
library(ggplot2)

read.cistrome.peak <- function(filepath, span = -1){
  # We need to remove rows with the extra chromosomes
  allowed_chr = c(paste("chr",1:22,sep=""),"chrX","chrY")
  df <- read.table(filepath)[ , c("V1", "V2", "V3", "V5")]
  colnames(df) <- c("chr","start","end","intensity")
  df <- filter(df, chr %in% allowed_chr)
  if(span != -1){
    # update start and end to use peak center, if span != -1
    peak_center <- (df$start + df$end) %/% 2 # integer division
    df$start <- peak_center - span
    df$end <- peak_center + span
  }
  setDT(df, key = names(df))
  setkey(df, "chr", "start","end")
  return(df)
}

get.bsite.cist.count <- function(peak_df, imads_df){
  sites_peak <- foverlaps(peak_df,imads_df)[, .(
    chr, peak.start = i.start, peak.end = i.end,
    intensity,
    pref
  )]

  sites_peak$count <- ifelse(is.na(sites_peak$pref), 0, 1)
  counted <- aggregate(sites_peak$count, by=list(sites_peak$chr, sites_peak$peak.start,
                                                 sites_peak$peak.end, sites_peak$intensity), FUN=sum)
  colnames(counted) <- c("chr","peak.start","peak.end","intensity","count")
  counted <- counted %>% arrange(chr,peak.start,peak.end)

  return(counted)
}

make.peaklen.dist.cistrome.plot <- function(lenvec1, outdir, chip_name=""){
  # hardcode upper limit to be 1500
  peaklen_hist_plot1 <- ggplot() +
    geom_histogram(aes(x = lenvec1), binwidth=50, colour = "white", fill = "cornflowerblue", size = 0.1) +
    labs(title="ChIP-seq peak length distribution",
         subtitle=chip_name,
         x="peak length") +
    coord_cartesian(xlim=c(min(lenvec1),1500)) +
    scale_x_continuous(breaks=seq(0, 1500, by = 100))
  ggsave(paste(outdir,"/",chip_name,"_peaklen_hist.pdf",sep=''), plot=peaklen_hist_plot1)

  peaklen_dist_plot <- ggplot() +
    geom_density(aes(lenvec1, fill="peaklen_dist"), alpha=0.5) +
    labs(title="ChIP-seq peak length distribution",
         subtitle=chip_name,
         x="peak length") +
    coord_cartesian(xlim=c(min(lenvec1),1500)) +
    scale_x_continuous(breaks=seq(0, 1500, by = 100))
    #theme(legend.position=c(.8,.8), legend.title=element_blank(), legend.spacing.x=unit(0.2,'cm')) # legend.title=element_blank()
  ggsave(paste(outdir,"/",chip_name,"_peaklen_dist.pdf",sep=''), plot=peaklen_dist_plot)
}

make.intensity.dist.plot <- function(count_df, outpath, chip_name="", maxsitesnum=8){
  # make boxplot, since we need pileup for the whole distribution, make a copy of the whole table
  # as group "all"
  merged_copy <- count_df
  merged_copy$count <- "all"
  merged_duplicated <- rbind(count_df,merged_copy)
  countvec <- c("all",as.character(sort(unique(count_df$count)))) # need this to order xlabel
  merged_duplicated$count <- factor(merged_duplicated$count, levels = countvec)

  pudist_plot <- ggplot(merged_duplicated, aes(x=count, y=intensity, fill=count)) +
    geom_boxplot() +
    #geom_jitter(width=0.1,alpha=0.2) +
    stat_summary(
      fun.data = count.n,
      geom = "text"
    ) +
    labs(title="Intensity distributions",
         subtitle=chip_name,
         y="intensity",
         x="#sites in a peak") +
    coord_cartesian(xlim=c(0,maxsitesnum+2)) +
    theme(legend.position = "none")
  ggsave(outpath)
}
