library(data.table)
library(plyr) 
library(dplyr) 
library(tidyr)
library(purrr)
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")

# kmer_alignment_file <- "runx1_kmer_alignment.txt"
# peak_file = '/Users/vincentiusmartin/Research/chip2gcPBM/resources/cistrome_files_runx1/44097.bed'
# corepos_pwm <- c(4,9)  # position of core in the pwm
# threshold <- 0.4

#' desc
#' @param peaklen_lists def
#' @keywords 
#' @export
#' @examples
#' run.kompas
run.kompas <- function(kmer_alignment_file, peak_file, corepos_pwm, threshold, 
                       zerobased=FALSE, genomever="hg19"){
  # table operation, we separate cells that have multiple values into different rows
  align_df <- fread(kmer_alignment_file, skip=7) %>%
    mutate(kposition = as.character(gsub("\\[|\\]", "",kposition))) %>%
    separate_rows(kposition, sep = ",") %>%
    mutate(kposition = as.numeric(kposition) + 1)
  # this is R so we use 1 indexing, that's why +1, TODO: don't make this assumption?
  #  align_df[align_df$kmer == "CAAGGGTC",]
  
  k <- nchar(align_df$kmer[1])
  core_len <- corepos_pwm[2] - corepos_pwm[1]
  center_pos <- k + 1
  
  # read the pwm file, only 4 rows
  pwm <- fread(kmer_alignment_file, skip=2, nrows=4)
  motif_start <- k + 1
  motif_len <- ncol(pwm) - 1
  motif_end <- k + motif_len
  
  #### checking here if needed 
  corepos_model <- corepos_pwm + k
  
  # reqkpos is the start region of every 8mer that encompasses the core
  if (k > core_len){
    reqkpos <- seq(corepos_model[2]-k+1, corepos_model[1])
  }else{
    reqkpos <- seq(corepos_model[1], corepos_model[2]-k)
  }
  
  # fill in the flank for threshold score reporting
  reqkpos <- seq(reqkpos[1]-1, reqkpos[length(reqkpos)]+1)
  align_df <- align_df %>% 
    filter(classified == 1 & kposition %in% reqkpos & Escore > threshold) %>%
    select(kmer, kposition, Escore)
  
  # read the peak files (e.g. encode, cistrome, etc)
  peak_df <- fread(peak_file, select = c(1:3))
  colnames(peak_df) <- c("chr", "start", "end")
  peak_df$start = peak_df$start + 1 # one indexing
  chr_order<-c(paste("chr",1:22,sep=""),"chrX","chrY")
  # remove non numeric chromosomes from the chromosomes column
  peak_df <- peak_df %>% filter(chr %in% chr_order)
  
  # figure out how long does a peak should be, to use them
  # vm: I think it's most likely that peak will fulfill this, so why do we need this
  min_peaklen <- ifelse(motif_len > k + 1, motif_len, k + 1)
  peak_df <- unique(peak_df[peak_df$end - peak_df$start > min_peaklen,]) # filter for short sequences the caller would have trouble with
  
  if (genomever == "hg19") {
    hs_gene <- BSgenome.Hsapiens.UCSC.hg19
  } else {
    hs_gene <- BSgenome.Hsapiens.UCSC.hg38
  }
  # getSeq(
  #   hs_gene,
  #   "chr17",
  #   start = 81285628,
  #   end = 81285852
  # )
  
  peak_sequences <- getSeq(
    hs_gene,
    peak_df$chr,
    start = peak_df$start,
    end = peak_df$end
  )
  
  # TODO: confirm indexing
  fwds <- as.character(peak_sequences, use.names=FALSE) # use_names = FALSE, no chr name
  rev_comps <- as.character(reverseComplement(peak_sequences), use.names=FALSE) 
  seq_lens <- width(peak_sequences)
  
  peak_list <- split(peak_df, seq(nrow(peak_df)))
  
  match_fwd <- lapply(fwds, kmer.match, escore_df = align_df, k = k)
  center_fwd <- mapply(find.center, match_fwd, orient="fwd", seq_lens, center_pos=center_pos, threshold=threshold)
  sites_fwd <- rbindlist(mapply(get.sites, center_fwd, peak_list, motif_len, orient='fwd', zerobased=zerobased, SIMPLIFY=FALSE)) %>%
    filter_all(all_vars(!is.na(.)))
  
  match_rc <- lapply(rev_comps, kmer.match, escore_df = align_df, k = k)
  center_rc <- mapply(find.center, match_rc, orient="rc", seq_lens, center_pos=center_pos, threshold=threshold)
  sites_rc <- rbindlist(mapply(get.sites, center_rc, peak_list, motif_len, orient='rc', zerobased=zerobased , SIMPLIFY=FALSE)) %>% 
    filter_all(all_vars(!is.na(.)))
  
  
  sites_all <- rbind(sites_fwd,sites_rc)
  sites_all$chr <- factor(sites_all$chr, levels=chr_order)
  sites_all <- sites_all %>% arrange(chr,start,end)
  
  return(sites_all)
}

#------------------


# https://stackoverflow.com/questions/35561641/find-all-possible-substrings-of-length-n
allsubstr <- function(x, n) substring(x, 1:(nchar(x) - n + 1), n:nchar(x))

#' def this is a hybrid of kmerMatch and findConsArrays in the original code
#' @param peaklen_lists def
#' @keywords 
#' @export
#' @examples
#' kmer.match()
kmer.match <- function(sequence, escore_df, k){
  kmers <- allsubstr(sequence, k)
  kmers_df <- data.frame("kmer" = kmers)
  scores <- join(x=kmers_df,y=escore_df,by="kmer",type = "left")
  
  #fill NA
  scores$kposition[is.na(scores$kposition)] <- 0
  scores$Escore[is.na(scores$Escore)] <- -0.5
  scores$pos <- rownames(scores)
  
  kpos <- with(rle(scores$Escore != -0.5), {
    ok <- values == TRUE & lengths > 1
    ends <- cumsum(lengths)
    starts <- ends - lengths + 1
    data.frame(starts, ends)[ok, ]
  })
  #kpos <- data.frame(start=c(2,41), end=c(5,45)) # testing purpose
  
  kresult <- map2_dfr(kpos$start, kpos$end, ~scores[.x:.y, ], .id = "group")
  # select(kmer, pos, kposition, orientation, Escore, group)
  return(kresult)
}

#' 
#' @param orient "fwd" or "rc"
#' @keywords 
#' @export
#' @examples
#' find.center()
find.center <- function(match_df, orient, seqlen, center_pos, threshold=0.4){
  if (nrow(match_df) == 0){
    return(NA)
  }
  match_df[c("kposition", "pos")] <- sapply(match_df[c("kposition", "pos")],as.numeric)
  md <- match_df %>% 
    group_by(group) %>% 
    arrange(kposition) %>%
    filter(row_number()==1) %>%
    mutate(
      center_site = center_pos - kposition + pos
    )
  if (orient == "rc"){
    md$center_site = seqlen - md$center_site + 1
  }
  return(sort(md$center_site))
}

get.sites <- function(center_list, peak_coord, motif_len, orient, zerobased=FALSE){
  #if(nrow) peak_coord
  seqlen <- peak_coord$end - peak_coord$start + 1
  centers <- unlist(center_list)
  start_pos <- if(orient == "fwd") peak_coord$start + centers - 1 else peak_coord$start + seqlen - centers - motif_len + 1
  ori <- if_else(orient == "fwd", "+", "-")
  sites <- data.frame("chr"= as.character(peak_coord$chr), "start"=start_pos, "end"=start_pos+motif_len-1,'ori'=ori)
  if(zerobased){
    sites$start <- sites$start - 1
  }
  return(sites)
}
