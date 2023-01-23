library(dplyr)
library(reshape2)
library(DESeq2)
library(SummarizedExperiment)
library(pheatmap)
rm(list=ls())

# script taken from /Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/RNAseq/Pan cancer/PBRM1_TCGA.R


#download data from GDC portal
dir <- "/Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/RNAseq/Pan cancer"
tcga_rnaseq_counts <- "/Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/RNAseq/Pan cancer/gdc_download_20200511_042315.478407/"
metadata <- read.table(paste0(dir,"/gdc_sample_sheet.KIRC.txt"),header = T,sep="\t")

setwd(tcga_rnaseq_counts)

file_list <- list.files(pattern="htseq.counts")

for (file in file_list){

  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=F, sep="\t")
  }

  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <- read.table(file, header=F, sep="\t")
    dataset <- cbind(dataset, file=temp_dataset[,2])
    rm(temp_dataset)
  }
}


rownames(dataset) <- dataset[,1]
dataset1 <- dataset[,-(1:2)]
dataset1 <- dataset1[-(60484:60488),]
colnames(dataset1) <- paste(file_list,".gz",sep="")
newnames <- metadata$Sample.ID[match(t(colnames(dataset1)),metadata$File.Name)]
colnames(dataset1) <- newnames


# PBRM1/VHL mutation data and cancer status data
# mutation assignment
setwd(dir)
broad <- read.table("broad.mit.edu__IlluminaGA_automated_DNA_sequencing_level2.maf",header=T, sep="\t", stringsAsFactors = F)
hgsc <- read.table("hgsc.bcm.edu__Mixed_DNA_Sequencing_level2.maf", header=T, sep="\t",quote = "", stringsAsFactors = F)
KIRC_ID <- unique(c(broad$Tumor_Sample_Barcode,hgsc$Tumor_Sample_Barcode))


Mutation <- array('NA',dim=c(ncol(dataset1),2))
colnames(Mutation) <- c('VHL','PBRM1')
Mutation[match(substr(KIRC_ID,1,16), colnames(dataset1)),]="WT"

broad_VHL <- broad[which(broad$Hugo_Symbol=="VHL"),]
broad_PBRM1 <- broad[which(broad$Hugo_Symbol == "PBRM1"),]
hgsc_VHL <- hgsc[which(hgsc$Hugo_Symbol=="VHL"),]
hgsc_PBRM1 <- hgsc[which(hgsc$Hugo_Symbol == "PBRM1"),]

Mutation[match(substr(hgsc_VHL$Tumor_Sample_Barcode, 1, 16),colnames(dataset1)), "VHL"]="MUT"
Mutation[match(substr(broad_VHL$Tumor_Sample_Barcode, 1, 16),colnames(dataset1)), "VHL"]="MUT"
Mutation[match(substr(hgsc_PBRM1$Tumor_Sample_Barcode, 1, 16),colnames(dataset1)), "PBRM1"]="MUT"
Mutation[match(substr(broad_PBRM1$Tumor_Sample_Barcode, 1, 16),colnames(dataset1)), "PBRM1"]="MUT"
rownames(Mutation) <- colnames(dataset1)


Sample.Type <- as.character(metadata$Sample.Type[match(rownames(Mutation),substr(metadata$Sample.ID,1,16))])
Mutation <- cbind(Mutation, Sample.Type)


# gene annotation
E2G <- read.delim2(file = "gencode.v22.annotation.gene.probeMap")
symbol <-  E2G$gene[match(rownames(dataset1),E2G$id)]
rowData <- data.frame(ID=rownames(dataset1), symbol = symbol)
rownames(rowData) <- rowData$ID

dataset1 <- data.frame(dataset1)
Mutation <- data.frame(Mutation)
rowData <- data.frame(rowData)

# create normalized count
dds <- DESeqDataSetFromMatrix(
  countData = dataset1,
  colData = Mutation,
  rowData = rowData,
  design = ~ Sample.Type )


dds <- DESeq(dds, )
dds <- dds[rowMeans(counts(dds)) > 1, ]
normcounts <- counts(dds, normalized = TRUE)
assays(dds)[["logcounts"]] <- log2(normcounts+1)

res <- results(dds, contrast=c("Sample.Type","Primary Tumor","Solid Tissue Normal"))
rowData(dds) <- cbind(rowData(dds), res)

current_dir <- "/Users/abc/Desktop/ASTAR/Projects/dCas9P300/screen results/"
setwd(current_dir)
saveRDS(dds,"data/TCGA.dds.rds")

dds <- readRDS(paste0(dir,"/TCGA.dds.rds"))


#Import 84 overlapping genes
gene_set <- read.delim2(file="data/New overlapped gene list.txt", stringsAsFactors=F, header=F, col.names=F)
gene_set <- gene_set[,1]
gene_set <- trimws(gene_set)

index <- as.numeric(na.omit(match(gene_set, rowData(dds)$symbol)))

dataset_set <- dds[index,]
dataset_set_counts <- assays(dataset_set)$logcounts

scale_counts <- t(scale(t(assays(dataset_set)$logcounts),center = T, scale = T))
annotation_col <-  data.frame(status=as.factor(dataset_set$Sample.Type), 
                              VHL=as.factor(dataset_set$VHL),
                              PBRM1=as.factor(dataset_set$PBRM1))
rownames(annotation_col) <- colnames(dataset_set)

rownames(scale_counts) <- rowData(dds)$symbol[index]

pdf("figures/TCGA.pdf", height = 12)  
pheatmap(scale_counts,
         show_colnames=FALSE,
         show_rownames=TRUE, 
         annotation_col=annotation_col, 
         breaks = seq(-2, 2, by = 4/31), 
         color = viridis::plasma(30))
dev.off()

# filter out genes overexpressed in T

gene_set <- intersect(rowData(dds)$symbol[rowData(dds)$log2FoldChange > 0.5], gene_set)

index <- as.numeric(na.omit(match(gene_set, rowData(dds)$symbol)))

dataset_set <- dds[index,]
dataset_set_counts <- assays(dataset_set)$logcounts

scale_counts <- t(scale(t(assays(dataset_set)$logcounts),center = T, scale = T))
annotation_col <-  data.frame(status=as.factor(dataset_set$Sample.Type), 
                              VHL=as.factor(dataset_set$VHL),
                              PBRM1=as.factor(dataset_set$PBRM1))
rownames(annotation_col) <- colnames(dataset_set)

rownames(scale_counts) <- rowData(dds)$symbol[index]

pdf("figures/TCGA.T.pdf", height = 9)  
pheatmap(scale_counts,
         show_colnames=FALSE,
         show_rownames=TRUE, 
         annotation_col=annotation_col, 
         breaks = seq(-2, 2, by = 4/31), 
         color = viridis::plasma(30))
dev.off()

