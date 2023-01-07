longrange <- read.table("../HiChIP/FitHiChIP/UNI3552_10kb.interactions_FitHiC_Q0.1_WashU.bed")

longrangeToInteract <- function(longrange, filename){
  source <- longrange[,1:3]
  target <- data.frame(do.call(rbind, lapply(longrange[,4], function(x){unlist(strsplit(x, split = ":|-|,"))})))
  span <- matrix(".", nrow=nrow(source), ncol=8)
  for (i in seq_len(nrow(source))){
    if (source[i,1] == target[i,1]){
      span[i,1] <- source[i,1]
      span[i,2] <- source[i,2]
      span[i,3] <- target[i,3]
    } else
      span[i, 1:3] <- source[i, 1:3]
  } 
  # score 
  span[,6] <- target[,4]
  span[,5] <- 0
  span[,8] <- 0
  
  final <- cbind(span, source, ".",".", target[,1:3],".",".")
  
  if (!is.null(filename)){
 
    header1 <- 'track type=interact name="interaction" description="An interact file" interactDirectional=true maxHeightPixels=200:100:50 visibility=full\n'
    header2 <- "#chrom  chromStart  chromEnd  name  score  value  exp  color  sourceChrom  sourceStart  sourceEnd  sourceName  sourceStrand  targetChrom  targetStart  targetEnd  targetName  targetStrand\n"
    cat(header1, file=filename)
    cat(header2, file=filename, append=TRUE)
    write.table(final, filename, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append=TRUE)
  }
 
  final
}


interact <- longrangeToInteract(longrange, filename = "../HiChIP/FitHiChIP/UNI3552_10kb.interactions_FitHiC_Q0.1.interact.txt")
