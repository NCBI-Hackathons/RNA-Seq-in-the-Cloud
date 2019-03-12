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
out_path <- "/home/jmcgirr/output/"
#
# Set project names
proj1 <- "ERP003613"
proj2 <- "ERP000546"

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

#####
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

res_ordered <- res[order(res$padj),]
res_ordered <- data.frame(res_ordered)
res_ordered$geneID <- rownames(res_ordered)
res_ordered$attribute_comparison <- paste(group1s[i],group2s[i],sep="_v_")
res_ordered$geneID <- rownames(res_ordered)

norm_cts <- data.frame(counts(dds, normalized=TRUE))
norm_cts$geneID <- rownames(norm_cts)
final_results <- merge(res_ordered,norm_cts, by = c("geneID"))

write.table(final_results,paste(group1s[i],"_vs_",group2s[i],"results.txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
}


