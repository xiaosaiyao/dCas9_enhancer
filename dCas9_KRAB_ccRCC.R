library("DESeq2")
library(pheatmap)

pcutoff <- 0.05
logFCcutoff <- -1



batch <- read.table('batch2.txt',header=TRUE, sep='\t',row.names=1,stringsAsFactors=FALSE)

# Perform differential analysis for A498
A498_matrix <- read.table('A498-krab-gained-enh.txt',header=TRUE, sep='\t',row.names=1,stringsAsFactors=FALSE)

A498_dds <- DESeqDataSetFromMatrix(countData = A498_matrix,
                              colData = batch,
                              design = ~ condition)


A498_dds <- DESeq(A498_dds, fitType="local")
A498_res <- results(A498_dds, contrast=c("condition","krab","Ctrl"))


# Plot MA plot
# In order to shade by pvalue, replace padj by pvalue just for plotting
A498_res.MA <- A498_res
A498_res.MA$padj <- A498_res.MA$pvalue
pdf("A498.R2.ma.pdf")
plotMA(A498_res.MA, 
       alpha = 0.05 ,
       ylim=c(-10,10), 
       main = "A498", 
       colSig = "red")
dev.off()


# Plot PCA to check
A498_vsd <- vst(A498_dds, blind=FALSE)
plotPCA(A498_vsd, intgroup=c("condition"))

write.table(A498_res,paste('Deseq2_A498_krab_R2.local.txt',sep=''),row.names=T, col.names=T,sep="\t")


# Plot heatmp 
A498_rld <- rlog(A498_dds,fitType='local')
A498_rld.sig <- assay(A498_rld)[which(A498_res$pvalue < pcutoff & A498_res$log2FoldChange < logFCcutoff), ]
A498_rld.sig <- data.frame(A498_rld.sig)
A498_res$enhancerID <- rownames(A498_res)

pdf("A498.heatmap.pdf", height=10, width = 3)
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

####################################################################################
# 786O

# Perform differential analysis for 786O

O786_matrix <- read.table('786O-krab-gained-enh.txt',header=T, sep='\t',row.names=1,stringsAsFactors=F)

O786_dds <- DESeqDataSetFromMatrix(countData = O786_matrix,
                              colData = batch,
                              design = ~ condition)




O786_dds <- DESeq(O786_dds, fitType="local")
O786_res <- results(O786_dds, contrast=c("condition","krab","Ctrl"),cooksCutoff=FALSE)

length(which(O786_combo$pvalue <0.05)) #1497/12330
length(which(O786_combo$padj <0.1)) #244/12330

# plot MA plot 
# In order to shade by pvalue, replace padj by pvalue just for plotting
O786_res.MA <- O786_res
O786_res.MA$padj <- O786_res.MA$pvalue
pdf("786O.ma.pdf")
plotMA(O786_res.MA, 
       alpha=0.05, 
       main="786O", 
       ylim=c(-10,10), 
       colSig = "red")
dev.off()

# Plot PCA to check
O786_vsd <- vst(O786_dds, blind=FALSE)
plotPCA(O786_vsd, intgroup=c("condition"))

# Plot cook's distance to check for outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(O786_dds)[["cooks"]]), range=0, las=2)

# Plot dispersion
plotDispEsts(O786_dds)
write.table(O786_res,paste('Deseq2_786O_krab_local.txt',sep=''),row.names=TRUE, col.names=TRUE, sep="\t")


# Plot heatmp 
O786_rld <- rlog(O786_dds, fitType='local')
O786_rld.sig <- assay(O786_rld)[which(O786_res$pvalue < pcutoff & O786_res$log2FoldChange < logFCcutoff), ]
O786_rld.sig <- data.frame(O786_rld.sig)
O786_res$enhancerID <- rownames(O786_res)


pdf("786O.heatmap.pdf", height=10, width = 3)
pheatmap(na.omit(O786_rld.sig),
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE, 
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

write.table(both_res,paste('Deseq2_786O_A498_krab_local.txt',sep=''),row.names=TRUE, col.names=TRUE,sep="\t")





# parse rownames to regions
gene_assign <- read.delim("great_gene_assign.txt")
gene1 <- unlist(lapply(strsplit(gene_assign$genes, split = ","), "[", 1))
gene1 <- unlist(lapply(strsplit(gene1, split = " "), "[", 1))

gene2 <- trimws(unlist(lapply(strsplit(gene_assign$genes, split = ","), "[", 2)))
gene2 <- unlist(lapply(strsplit(gene2, split = " "), "[", 1))
gene_assign$gene1 <- gene1
gene_assign$gene2 <- gene2

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

sig_genes <- sort(intersect(unique(c(O786_rld.sig$gene1, O786_rld.sig$gene2)), 
        unique(c(A498_rld.sig$gene1, A498_rld.sig$gene2))))



sort(table(c(O786_rld.sig$gene1, O786_rld.sig$gene2, A498_rld.sig$gene1, A498_rld.sig$gene2)))

test <- read.delim("overlap.test.txt")

intersect(unlist(test), sig_genes)