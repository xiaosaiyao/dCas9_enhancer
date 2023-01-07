# Differential analysis of dCas9_KRAB + gRNA targeting top 500 gained enhancer regions in A498 and 786O
# Top enhancers based on Yao, Cancer Discovery 2017
# Experiment performed by Tyler Klann of Charles Gersbach lab of Duke University

library(DESeq2)
library(pheatmap)

# Define differential analysis pvalue and log fold change cutoff
# Note that pvalue is used in place of p-adjusted to allow more exclusive secondary screening later on
pcutoff <- 0.05
logFCcutoff <- -1


batch <- read.table('data/batch2.txt', header=TRUE, sep='\t',row.names=1,stringsAsFactors=FALSE)

# Perform differential analysis for A498
A498_matrix <- read.table('data/A498-krab-gained-enh.txt',
                          header=TRUE, 
                          sep='\t',
                          row.names=1,
                          stringsAsFactors=FALSE)

A498_dds <- DESeqDataSetFromMatrix(countData = A498_matrix,
                                   colData = batch,
                                   design = ~ condition)


A498_dds <- DESeq(A498_dds, fitType="local")
A498_res <- results(A498_dds, contrast=c("condition","krab","Ctrl"))

length(which(A498_res$pvalue <0.05 & A498_res$log2FoldChange <0)) #479/12330
length(which(A498_res$pvalue <0.05 & A498_res$log2FoldChange >0)) #517/12330
length(which(A498_res$padj <0.1)) # 53/12330

# Plot MA plot
# In order to shade significant regions by pvalue, replace padj by pvalue just for plotting
A498_res.MA <- A498_res
A498_res.MA$padj <- A498_res.MA$pvalue
pdf("figures/A498.R2.ma.pdf")
plotMA(A498_res.MA, 
       alpha = 0.05 ,
       ylim=c(-5,5), 
       main = "A498", 
       colSig = "red")
dev.off()


# Plot PCA to check
A498_vsd <- vst(A498_dds, blind=FALSE)

pdf("figures/PCA.A498.pdf") # dCas9_KRAB replicates not tightly clustered
plotPCA(A498_vsd, intgroup=c("condition"))
dev.off()



# Plot cook's distance to check for outliers
pdf("figures/cookdistance.A498.pdf")
par(mar=c(8,5,2,2))
boxplot(log10(assays(A498_dds)[["cooks"]]), range=0, las=2)
dev.off()

# Plot dispersion
pdf("figures/dispersion.A498.pdf")
plotDispEsts(A498_dds)
dev.off()

# Plot heatmap 
A498_rld <- rlog(A498_dds,fitType='local')
A498_rld.sig <- assay(A498_rld)[which(A498_res$pvalue < pcutoff & A498_res$log2FoldChange < logFCcutoff), ]
A498_rld.sig <- data.frame(A498_rld.sig)
A498_res$enhancerID <- rownames(A498_res)

pdf("figures/A498.heatmap.pdf", height=5, width = 3)
pheatmap(na.omit(A498_rld.sig),
         scale = "row",
         show_colnames = TRUE,
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE,  
         breaks = seq(-2, 2, by = 4/31), 
         color = viridis::plasma(30),
         border_color = "NA")
dev.off()

write.table(A498_res,paste('results/Deseq2_A498_krab_R2.local.txt',sep=''),row.names=T, col.names=T,sep="\t")


####################################################################################
# 786O

# Perform differential analysis for 786O

O786_matrix <- read.table('data/786O-krab-gained-enh.txt',header=T, sep='\t',row.names=1,stringsAsFactors=F)

O786_dds <- DESeqDataSetFromMatrix(countData = O786_matrix,
                              colData = batch,
                              design = ~ condition)




O786_dds <- DESeq(O786_dds, fitType="local")
O786_res <- results(O786_dds, contrast=c("condition","krab","Ctrl"),cooksCutoff=FALSE)


length(which(O786_res$pvalue <0.05 & O786_res$log2FoldChange <0)) #471/12330
length(which(O786_res$pvalue <0.05 & O786_res$log2FoldChange >0)) #1026/12330

length(which(O786_res$padj <0.1)) #244/12330

# plot MA plot 
# In order to shade by pvalue, replace padj by pvalue just for plotting
O786_res.MA <- O786_res
O786_res.MA$padj <- O786_res.MA$pvalue
pdf("figures/786O.ma.pdf")
plotMA(O786_res.MA, 
       alpha=0.05, 
       main="786O", 
       ylim=c(-5,5), 
       colSig = "red")
dev.off()

# Plot PCA to check
pdf("figures/PCA.786O.pdf") # dCas9_KRAB replicates not tightly clustered
O786_vsd <- vst(O786_dds, blind=FALSE)
plotPCA(O786_vsd, intgroup=c("condition"))
dev.off()

# Plot cook's distance to check for outliers
pdf("figures/cookdistance.786O.pdf")
par(mar=c(8,5,2,2))
boxplot(log10(assays(O786_dds)[["cooks"]]), range=0, las=2)
dev.off()

# Plot dispersion
pdf("figures/dispersion.786O.pdf")
plotDispEsts(O786_dds)
dev.off()

write.table(O786_res,paste('results/Deseq2_786O_krab_local.txt',sep=''),row.names=TRUE, col.names=TRUE, sep="\t")


# Plot heatmap 
O786_rld <- rlog(O786_dds, fitType='local')
O786_rld.sig <- assay(O786_rld)[which(O786_res$pvalue < pcutoff & O786_res$log2FoldChange < logFCcutoff), ]
O786_rld.sig <- data.frame(O786_rld.sig)
O786_res$enhancerID <- rownames(O786_res)


pdf("figures/786O.heatmap.pdf", height=5, width = 3)
pheatmap(na.omit(O786_rld.sig),
         scale = "row",
         show_colnames = TRUE,
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE,  
         breaks = seq(-2, 2, by = 4/31), 
         color = viridis::plasma(30),
         border_color = "NA")
dev.off()

#########################################################
# Combine results
colnames(O786_res) <- paste0("786O_",colnames(O786_res))
colnames(O786_matrix) <- paste0("786O_",colnames(O786_matrix))
colnames(A498_res) <- paste0("A498_",colnames(A498_res))
colnames(A498_matrix) <- paste0("A498_",colnames(A498_matrix))
both_res <- cbind(O786_res, O786_matrix[rownames(O786_res),], A498_res,  A498_matrix[rownames( A498_res),])
both_res <- data.frame(both_res)
head(both_res)

write.table(both_res,paste('results/Deseq2_786O_A498_krab_local.txt',sep=''),row.names=TRUE, col.names=TRUE,sep="\t")





# parse rownames to regions
gene_assign <- read.delim("results/great_gene_assign.txt")
gene1 <- unlist(lapply(strsplit(gene_assign$genes, split = ","), "[", 1))
gene1 <- unlist(lapply(strsplit(gene1, split = " "), "[", 1))

gene2 <- trimws(unlist(lapply(strsplit(gene_assign$genes, split = ","), "[", 2)))
gene2 <- unlist(lapply(strsplit(gene2, split = " "), "[", 1))
gene_assign$gene1 <- gene1
gene_assign$gene2 <- gene2

# Merge additional information in significant regions
A498_rld.sig$enhancerID <- rownames(A498_rld.sig)
A498_rld.sig <- merge(A498_rld.sig, gene_assign)
A498_res <- data.frame(A498_res)
A498_res$enhancerID <- rownames(A498_res)
A498_rld.sig <- merge(A498_rld.sig, A498_res)



O786_rld.sig$enhancerID <- rownames(O786_rld.sig)
O786_rld.sig <- merge(O786_rld.sig, gene_assign)
O786_res <- data.frame(O786_res)
O786_res$enhancerID <- rownames(O786_res)
O786_rld.sig <- merge(O786_rld.sig, O786_res)

#export significant regions as bed files
A498.sig.bed <- do.call(rbind, lapply(A498_rld.sig$enhancerID, function(x){
        unlist(strsplit(x, split = "_"))}))
write.table(A498.sig.bed, "results/A498.sig.bed", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(A498_rld.sig, "results/A498_rld.sig.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(O786_rld.sig, "results/O786_rld.sig.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

O786.sig.bed <- do.call(rbind, lapply(O786_rld.sig$enhancerID, function(x){
        unlist(strsplit(x, split = "_"))}))
write.table(O786.sig.bed, "results/O786.sig.bed", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Intersect regions that are significant in both A498 and 786O
# gRNA
overlap_euler_regions <- eulerr::euler(combinations = list(
        O786 = O786_rld.sig$enhancerID,
        A498 = A498_rld.sig$enhancerID))

pdf("figures/venn.regions.gRNA.overlap.pdf", height =3, width = 3)
plot(overlap_euler_regions, quantities = list(type = c( "counts")))
dev.off()

# enhancers
overlap_euler_enhancers <- eulerr::euler(combinations = list(
        O786 = unique(O786.sig.bed[,4]),
        A498 = unique(A498.sig.bed[,4])))

pdf("figures/venn.enhancers.overlap.pdf", height =3, width = 3)
plot(overlap_euler_enhancers, quantities = list(type = c( "counts")))
dev.off()

# Intersect genes that are significant in both A498 and 786O

O786_sig_genes <- unique(c(O786_rld.sig$gene1, O786_rld.sig$gene2))
A498_sig_genes <- unique(c(A498_rld.sig$gene1, A498_rld.sig$gene2))
sig_genes <- sort(intersect(O786_sig_genes, A498_sig_genes))

overlap_euler <- eulerr::euler(combinations = list(
        O786 = O786_sig_genes,
        A498 = A498_sig_genes))

pdf("figures/venn.gene.overlap.pdf", height =3, width = 3)
plot(overlap_euler, quantities = list(type = c( "counts")))
dev.off()
sort(table(c(O786_rld.sig$gene1, O786_rld.sig$gene2, A498_rld.sig$gene1, A498_rld.sig$gene2)))

# Load the genes that were tested by Dylan and perform intersect
test <- read.delim("results/overlap.test.txt", header = FALSE)
intersected.genes <- intersect(unlist(test), sig_genes)

write.table(intersected.genes, "results/sig.genes.tested.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)