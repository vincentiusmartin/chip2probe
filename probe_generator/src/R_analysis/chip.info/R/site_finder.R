
# https://stackoverflow.com/questions/56957519/how-to-get-the-index-of-a-string-of-certain-length-that-is-made-from-a-list-of-c/56958112#56958112
seq2i <- function(seq){
  #Convert string to a vector of letters
  strletters <- unlist(strsplit(seq, ""))

  facts <- factor(strletters, levels=c("A", "C", "G", "T"))
  nums <- as.integer(facts)-1

  #create list of multipliers
  multipliers <- 4**((length(nums)-1):0)
  intrep <- sum(multipliers*nums) + 1

  return(intrep)
}

i2seq <- function(int){
  seqlen <- ceiling(log(int-1, base = 4))
  divisors <- 4 **((seqlen-1):0)
  quotients <- (int-1) %/% divisors
  remainders <- quotients %% 4
  bases <- c('A','C','G','T')
  seq <- paste(bases[as.factor(remainders)],collapse='')
  return(seq)
}
