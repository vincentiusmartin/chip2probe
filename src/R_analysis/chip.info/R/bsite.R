library(data.table)
library("BSgenome.Hsapiens.UCSC.hg19")

imads_models <- c()

read.imads.bed <- function(bed_path){
  imads_sites <- read.csv(bed_path,sep="\t",header=FALSE)[,1:4]
  colnames(imads_sites) <- c("chr","start","end","pref")
  setDT(imads_sites, key = names(imads_sites))
  setkey(imads_sites, "chr", "start","end")
  return(imads_sites)
}

remove.nosites.peak <- function(nrwp_df, imads_df){
  # remove peak with no match
  peak_wsites <- foverlaps(imads_df, nrwp_df, type="any", nomatch=NULL) %>% select(-c(start,end,pref))
  peak_wsites <- unique(peak_wsites)
  setDT(peak_wsites, key = names(peak_wsites))
  return(peak_wsites)
}

get.bsite.count <- function(nrwp_df, imads_df, summit){
  #setDT(nrwp_df, key = names(nrwp_df))
  #setkey(nrwp_df, "chr", "start","end")
  sites_peak <- foverlaps(nrwp_df,imads_df)[, .(
    chr, peak.start = i.start, peak.end = i.end,
    summit,
    pref
  )]

  sites_peak$count <- ifelse(is.na(sites_peak$pref), 0, 1)
  counted <- aggregate(sites_peak$count, by=list(sites_peak$chr,sites_peak$peak.start,sites_peak$peak.end), FUN=sum)
  colnames(counted) <- c("chr","peak.start","peak.end","count")

  return(counted)
}

get.bsites.in.peak <- function(pileups_df, imads_df, score_colname = "pileups"){
  # TODO: need to think about score_colname
  sites_selected <- foverlaps(imads_df, pileups_df, nomatch=NULL, type="within")
  # order the chromosome
  chr_order<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  sites_selected$chr <- factor(sites_selected$chr, levels=chr_order)
  sites_selected <- sites_selected[order(sites_selected$chr),] %>%
    rename(start = "bsite.start", end = "bsite.end")

  sites_in_peak_df <- sites_selected %>%
    group_by(chr, peak.start, peak.end) %>%
    mutate(distance = (bsite.start + bsite.end) / 2 - lag((bsite.start + bsite.end)/2),
           bsite1.start = lag(bsite.start),
           bsite1.end = lag(bsite.end),
           bsite2.start = bsite.start,
           bsite2.end = bsite.end,
           bsite1.pref = lag(pref),
           bsite2.pref = pref) %>%
    filter(!is.na(distance)) %>%
    #select(-c(count, start, end, pref)) %>%
    select(c(chr, peak.start, peak.end, score_colname, distance,
             bsite1.start, bsite1.end, bsite2.start, bsite2.end,
             bsite1.pref, bsite2.pref)) %>%
    select(-distance,everything()) %>% # just to put distance as the last col
    mutate_at(score_colname, round, 4)

  return(sites_in_peak_df)
}

get.bsites.within.range <- function(sites_in_peak_df, min_dist, max_dist, probeseq_flank = 0){
  # TODO: filter chrM?
  # Get binding sites in min_dist .. max_dist
  # if no sites are within distance, a warning message would show up:
  # Factor `chr` contains implicit NA, consider using `forcats::fct_explicit_na`
  sites_within_range <- sites_in_peak_df %>%
    filter(distance >= min_bsite_dist & distance <= max_bsite_dist)  %>%
    rowwise() %>% # need to calculate each row separately
    mutate(
      seqstart = ((bsite1.start + bsite1.end) / 2) - ceiling((probe_size - distance) / 2) + 1,
      seqend =  ((bsite2.start + bsite2.end) / 2) + floor((probe_size - distance) / 2),
      sequence = toString(
        getSeq(
          Hsapiens,
          chr,
          start = seqstart,
          end = seqend
        )
      ),
      flank_left = toString(
        getSeq(
          Hsapiens,
          chr,
          start = seqstart - probeseq_flank,
          end = seqstart - 1
        )
      ),
      flank_right = toString(
        getSeq(
          Hsapiens,
          chr,
          start = seqend +1 ,
          end = seqend + probeseq_flank
        )
      )
    )
    # %>% select(-c(seqstart,seqend)) # made just as a tmp variable so remove before returning
  return(sites_within_range)
}

# ---------- Predict Binding Site ----------

load.imads.models <- function(model_paths){
  models <- lapply(model_paths,readRDS)
}

# Installing libsvm in R: https://www.csie.ntu.edu.tw/~cjlin/libsvm/

# predict.binding <- function(sequences){
#
# }

# ---------- E-score ----------

read.escore.file <- function(escore_short_path, escore_map_path){
  eshort_df <- read.csv(escore_short_path, header=F)$V1
  emap_df <- read.csv(escore_map_path, header=T)$idx_short
  elong_df <- eshort_df[emap_df]
  return(elong_df)
}

# this should be in main, but I couldn't find svm_load_model implementation, so let's move it to python
# escore_short_path <- "/Users/vincentiusmartin/Research/chip2gcPBM/escores/ets1_escores.txt"
# escore_map_path <- "/Users/vincentiusmartin/Research/chip2gcPBM/escores/index_short_to_long.csv"
#
# imads_model_paths <- c("/Users/vincentiusmartin/Research/chip2gcPBM/imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAA_1a2a3mer_format.model",
#                        "/Users/vincentiusmartin/Research/chip2gcPBM/imads_files/models/ets1/ETS1_100nM_Bound_filtered_normalized_transformed_20bp_GGAT_1a2a3mer_format.model")

#escore_df <- read.escore.file(escore_short_path, escore_map_path)
# eshort_df <- read.escore.file(escore_short_path,escore_map_path)
#
# seq <- "GCAGAAGGGGCATCCGGGCAGCTTCCTGGTGCACAT"
# substrs <- substring(seq, 1:(nchar(seq)-8+1), 8:nchar(seq))
# eidx <- unlist(lapply(substrs, seq2i))
# escores <- elong_df[eidx]
# seqvec <- unlist(strsplit(seq,split=""))
# xindex <- c(4:32) # for escore, 8/2 to 36-4
#
# g <- ggplot() +
#   # need to make x as factor to make it into a discrete scale, also need to group
#   geom_line(aes(y=escores,x=factor(xindex,levels=c(1:length(seqvec))),group=1),colour="darkorange2",size=1.5) +
#   # need drop to use factor levels in x axis
#   scale_x_discrete(breaks=seq_along(seqvec),labels=seqvec,drop=FALSE) +
#   labs(title="Test",
#        subtitle=chip_name,
#        y="escores",
#        x="sequence")
# plot(g)
