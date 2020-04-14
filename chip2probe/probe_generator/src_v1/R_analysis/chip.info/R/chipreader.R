library(data.table)
library(dplyr)

#' This function does something
#' @param
#' @return
#' @export
#' @examples
#'
read.pileup <- function(pileup_path){
  pu <- read.table(pileup_path)
  colnames(pu) <- c("chr","start","end","pileup")
  setDT(pu, key = names(pu))
  setkey(pu, "chr","start","end")
  return(pu)
}

read.narrow.peak <- function(nrwp_path, span = -1){
  nrwp <- read.table(nrwp_path)[ , c("V1", "V2", "V3", "V10")]
  colnames(nrwp) <- c("chr","start","end","summit")
  if(span != -1){
    nrwp$end = nrwp$start + nrwp$summit + span
    nrwp$start = nrwp$start + nrwp$summit - span
  }
  setDT(nrwp, key = names(nrwp))
  setkey(nrwp, "chr", "start","end")
  return(nrwp)
}

get.pileups <- function(nrwp_df, pu_df, logpileup=FALSE, join_type="any"){
  pu_nrwp <- foverlaps(pu_df, nrwp_df, nomatch=NULL, type=join_type)[, .(
    chr, peak.start = start, peak.end = end,
    pileup.start = i.start, pileup.end = i.end,
    pileup,
    summit
  )]

  # Handle edge cases here, consider partial reads by bordering with both ends of the peak
  pu_nrwp$pileup.start <- ifelse(pu_nrwp$pileup.start < pu_nrwp$peak.start, pu_nrwp$peak.start, pu_nrwp$pileup.start)
  pu_nrwp$pileup.end <- ifelse(pu_nrwp$pileup.end > pu_nrwp$peak.end, pu_nrwp$peak.end, pu_nrwp$pileup.end)

  agg <- pu_nrwp %>%
    group_by(chr,peak.start,peak.end) %>%
    summarise(pileups = sum((pileup.end - pileup.start) * pileup))

  if(logpileup){
    # +1 to avoid -inf on log
    agg$pileups = log(agg$pileups+1)
  }

  return(agg)
}

write.peak.info <- function(nrwp_preidr, nrwp_postidr, path){
  cat(c(paste("#peaks before idr", nrow(nrwp_preidr)),paste("#peaks after idr", nrow(nrwp_postidr))))
  filconn<-file(path)
  writeLines(c(paste("#peaks before idr", length(nrwp_preidr)),paste("#peaks after idr", length(nrwp_postidr))), filconn)
  close(filconn)
}
