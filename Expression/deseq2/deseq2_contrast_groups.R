#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list (

    make_option (c("-c","--counts"),
        default="/home/tommerschwarz/data/ERP000546_genecounts.txt",
        help="The gene counts file for each comparable run [default %default]"),

    make_option (c("-a","--attributes"),
        default="/home/tommerschwarz/data/ERP000546_attributes.txt",
        help="The attributes file for info on each run [default %default]"),

    make_option (c("-t","--type"),
        help="A string matching a sample type in the 'type' column of the attributes file [default %default]"),
    
    make_option (c("-ols","--outliers"),
                 default=FALSE,
                 action="store_true",
                 help="Use flag to include outliers. [default %default]"),
    
    make_option (c("-o","--outdir"),
                 default="/home/jmcgirr/output/",
                 help="The attributes file for info on each run [default %default]")
    )


opt  <- parse_args(
    OptionParser(#usage= "usage: %prog [options]",
        option_list=option_list)
        )
print(opt$outliers)
print(opt$type)
suppressPackageStartupMessages(library(DESeq2))

############################################
###### Set paths and load 2 dataframes #####
############################################

# cluster jmcgirr@aa-rnaseq-clustering
# Set paths for two dataframes created by Tommer
# 1. read counts for each comparable run
counts_path <- opt$counts
#
# 2. attributes for each run
atts_path <- opt$attributes
#
# Set output path for PCA plots, log fold-change tables, normalized counts
out_path <- opt$outdir

show.warnings=FALSE
dir.create(out_path)
show.warnings=TRUE

setwd(out_path)

# local
#counts_path <- "C:/Users/jmcgirr/Desktop/joint_genecounts.txt"
#atts_path <- "C:/Users/jmcgirr/Desktop/joint_metadata.txt"

# What groups do we want to compare?
# Which column in the attributes should be compared?
att_cols <- c("condition")
group1s <- c("control")
group2s <- c("case")

#####
############################################
### Create DESeq2 object and run contrast ##
############################################

atts <- read.table(atts_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
atts <- atts[which(atts$type == opt$type & atts$is_outlier == "False"),]
counts <- read.table(counts_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
rownames(counts) <- counts[,1]
counts <- counts[,-1]
keeps <- names(counts)[(names(counts) %in% atts$Run)]
counts <- counts[, keeps]
rownames(atts) <- atts[,1]
atts <- atts[,-1]
atts <- as.matrix(atts)

if (opt$outliers)
{
  atts <- read.table(atts_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  counts <- read.table(counts_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  atts <- atts[which(atts$type == opt$type),]
  rownames(counts) <- counts[,1]
  counts <- counts[,-1]
  keeps <- names(counts)[(names(counts) %in% atts$Run)]
  counts <- counts[, keeps]
  rownames(atts) <- atts[,1]
  atts <- atts[,-1]
  atts <- as.matrix(atts)
}

ncol(counts)
nrow(atts)

# create DESeq Object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = atts,
                              design= ~condition+project)
dds$condition <- factor(dds$condition, levels = c("control","case"))

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 2 ) >= 2
dds <- dds[idx,]
dds <- DESeq(dds)
rld <- vst(dds)

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

write.table(final_results,paste(group1s[i],"_vs_",group2s[i],"_results.txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
}


