library(Repitools)
library(rtracklayer)

# split hg19 into equal bins of 100KB
hg19_size <- read.table("data/hg19.chrom.sizes")
hg19_v <- hg19_size[,2]
names(hg19_v) <- hg19_size[,1]
hg19_blocks <- genomeBlocks(hg19_v, chrs = names(hg19_v), width = 100000)

# load significant depleted gRNA regions
A498.sig <- read.delim("results/A498.sig.bed", header = FALSE, sep = " ")
colnames(A498.sig) <- c("chr","start","end","enhancerID","gRNA")
A498.sig <- makeGRangesFromDataFrame(A498.sig, keep.extra.columns=TRUE)

O786.sig <- read.delim("results/O786.sig.bed", header = FALSE, sep = " ")
colnames(O786.sig) <- c("chr","start","end","enhancerID","gRNA")
O786.sig <- makeGRangesFromDataFrame(O786.sig, keep.extra.columns=TRUE)

both.sig <- c(O786.sig, A498.sig)
both.sig <- reduce(both.sig)

# count depleted gRNAs in each bin 
genomic_counts <- findOverlaps(both.sig, hg19_blocks)

bin_counts <- sort(table(as.data.frame(genomic_counts)[,2]),decreasing = TRUE)
sig.blocks <- hg19_blocks[as.numeric(names(bin_counts))]
sig.blocks$counts <- bin_counts
rtracklayer::export(object = sig.blocks, con = "results/sig.blocks.bed", format = "BED")
