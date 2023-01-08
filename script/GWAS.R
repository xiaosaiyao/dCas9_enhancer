library(rtracklayer)

# Load all RCC GWAS snps
GWAS <- read.csv( "data/efotraits_EFO_0000681-associations-2023-01-7.hg19.csv")
GWAS <- makeGRangesFromDataFrame(GWAS, keep.extra.columns=TRUE)

# Load all enhancers with depleted gRNAs
sig_enhancers <- read.table("results/sig.enhancers.A498.786O.txt")
colnames(sig_enhancers) <- c("chr","start","end", "name")
sig_enhancers <- makeGRangesFromDataFrame(sig_enhancers, keep.extra.columns=TRUE)


GWAS_overlap <- findOverlaps(GWAS, sig_enhancers, maxgap = 50000)
export(GWAS[GWAS_overlap@from], "results/GWAS.overlap.bed", "BED")
export(sig_enhancers[GWAS_overlap@to], "results/GWAS.enhancers.overlap.bed", "BED")