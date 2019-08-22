library(ggplot2)

# theme_set(theme_bw())

make.peaklen.dist.plot <- function(lenvec1, lenvec2, chip_name="", dir){
  # hardcode upper limit to be 1500
  peaklen_hist_plot1 <- ggplot() +
    geom_histogram(aes(x = lenvec1), binwidth=50, colour = "white", fill = "cornflowerblue", size = 0.1) +
    labs(title="ChIP-seq peak length r1",
         subtitle=chip_name,
         x="peak length") +
    coord_cartesian(xlim=c(0,1500)) +
    scale_x_continuous(breaks=seq(0, 1500, by = 100))
  ggsave(paste(dir,"/peaklen_hist_before_idr.pdf",sep=''), plot=peaklen_hist_plot1)

  peaklen_hist_plot2 <- ggplot() +
    geom_histogram(aes(x = lenvec2), binwidth=50, colour = "white", fill = "cornflowerblue", size = 0.1) +
    labs(title="ChIP-seq peak length r2",
         subtitle=chip_name,
         x="peak length") +
    coord_cartesian(xlim=c(0,1500)) +
    scale_x_continuous(breaks=seq(0, 1500, by = 100))
  ggsave(paste(dir,"/peaklen_hist_after_idr.pdf",sep=''), plot=peaklen_hist_plot2)

  peaklen_dist_plot <- ggplot() +
    geom_density(aes(lenvec1, fill="before idr"), alpha=0.5) +
    geom_density(aes(lenvec2, fill="after idr"),alpha=0.5) +
    labs(title="ChIP-seq peak length distribution",
         subtitle=chip_name,
         x="peak length") +
    coord_cartesian(xlim=c(min(lenvec1),1500)) +
    scale_x_continuous(breaks=seq(0, 1500, by = 100)) +
    #labs(fill="Peak length distribution") +
    theme(legend.position=c(.8,.8), legend.title=element_blank(), legend.spacing.x=unit(0.2,'cm')) # legend.title=element_blank()
  #peaklen_dist_plot
  ggsave(paste(dir,"/peaklen_dist.pdf",sep=''))
}

make.corr.plot <- function(xvals, yvals, path, chip_name = ""){
  # pdf(path)
  # plot(merged$pileups.x, merged$pileups.y, cex=0.2, col = "blue")
  # abline(coef = c(0,1), col = "red")
  # pu_cor <- cor(x, y)
  # title(paste("Correlation",title))
  # legend("topleft",legend=sprintf("R² = %.2f",pu_cor^2), bty ="n", pch=NA)
  # dev.off()

  # use ggplot instead

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

count.n <- function(x) {
  return(c(y = median(x)*1.03, label = length(x)))
}

make.pileup.dist.plot <- function(bsites_ctpile, outpath, chip_name){
  # make boxplot, since we need pileup for the whole distribution, make a copy of the whole table
  # as group "all"
  merged_copy <- bsites_ctpile
  merged_copy$count <- "all"
  merged_duplicated <- rbind(bsites_ctpile,merged_copy)
  countvec <- c("all",as.character(sort(unique(bsites_ctpile$count)))) # need this to order xlabel
  merged_duplicated$count <- factor(merged_duplicated$count, levels = countvec)
  pudist_plot <- ggplot(merged_duplicated, aes(x=count, y=pileups, fill=count)) +
    geom_boxplot() +
    #geom_jitter(width=0.1,alpha=0.2) +
    stat_summary(
      fun.data = count.n,
      geom = "text"
    ) +
    labs(title="Pileup distributions",
         subtitle=chip_name,
         y="log(pileup)",
         x="#sites in a peak") +
    theme(legend.position = "none")
  ggsave(outpath)

  # count here
  #count_df <- as.data.frame(table(merged_duplicated$count))
  #write.table(count_df,file=paste(outpath,"/bsite_count_",span,".tsv",sep=''),sep="\t",row.names = FALSE, quote = FALSE)
}
