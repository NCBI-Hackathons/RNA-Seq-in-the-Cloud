library(DESeq2)

############################################
###### Set paths and load 2 dataframes #####
############################################

# cluster jmcgirr@aa-rnaseq-clustering
# Set paths for two dataframes created by Tommer
# 1. read counts for each comparable run
counts_path <- "/home/tommerschwarz/data/ERP000546_genecounts.txt"
#
# 2. attributes for each run
atts_path <- "/home/tommerschwarz/data/ERP000546_attributes.txt"
#
# Set output path for PCA plots, log fold-change tables, normalized counts

# local
#counts_path <- "C:/Users/jmcgirr/Documents/GitHub/RNA-Seq-in-the-Cloud/Expression/deseq2/data/ERP000546_genecounts.txt"
#atts_path <- "C:/Users/jmcgirr/Documents/GitHub/RNA-Seq-in-the-Cloud/Expression/deseq2/data/ERP000546_attributes.txt"

# Set column in attribute(s) file to include in DESeq objects
design <- "type"

# What groups do we want to compare?
# Which column in the attributes should be compared?
att_cols <- c("type")
group1s <- c("a")
group2s <- c("b")

############################################
###### Create DESeq2 object and run ########
############################################

counts <- read.table(counts_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
rownames(counts) <- counts[,1]
counts <- counts[,-1]
atts <- read.table(atts_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
atts <-as.matrix(read.table(atts_path ,header = TRUE,row.names=1))

ncol(counts)
nrow(atts)

# create DESeq Object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = atts,
                              design= ~type)

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 2 ) >= 2
dds <- dds[idx,]
dds <- DESeq(dds)
rld <- vst(dds)

i <- 1
for (i in c(1:length(att_cols)))
{
att_cols[i]
group1s[i]
group2s[i]
  
# contrast groups
res <- results(dds, contrast=c(att_cols[i],group1s[i],group2s[i]))
#results(dds, contrast=c(att_cols[i]))

pca_plot <- plotPCA(rld, intgroup=c(att_cols[i]))

#tiff(paste(out_path,"pca.tiff",sep = ""), width = 6, height = 6, units = 'in', res = 500)
print(pca_plot)
#dev.off()

resOrdered <- res[order(res$padj),]
res_ordered <- as.data.frame(resOrdered)
res_ordered$geneID <- rownames(res_ordered)
print(head(resOrdered))

#resLFC <- lfcShrink(dds, coef=2, res=res)

total_genes <- nrow(res_ordered)
de_total <- nrow(res_ordered[which(res_ordered$padj < 0.05),])
de_up <- (nrow(res_ordered[which(res_ordered$log2FoldChange > 0 & res_ordered$padj < 0.05),]))/total_genes
de_dn <- (nrow(res_ordered[which(res_ordered$log2FoldChange < 0 & res_ordered$padj < 0.05),]))/total_genes
prop_de <- de_total/total_genes
total_genes_plot <- paste(total_genes, "genes", sep = " ")

ma_plot <- plotMA(res, colSig = "darkred", colNonSig = "grey")
#tiff(paste(out_path,"de_plot.tiff"sep = ""), width = 3.5, height = 3.5, units = 'in', res = 500)
print(ma_plot)
#dev.off()

}