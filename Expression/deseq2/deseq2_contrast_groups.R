

suppressPackageStartupMessages(library(optparse))

option_list <- list (

    make_option (c("-c","--counts"),
        default="/home/tommerschwarz/data/ERP000546_genecounts.txt",
        help="The gene counts file for each comparable run [default %default]"),

    make_option (c("-a","--attributes"),
        default="/home/tommerschwarz/data/ERP000546_attributes.txt",
        help="The attributes file for info on each run [default %default]"),

    make_option (c("-o","--outdir"),
        default="/home/jmcgirr/output/",
        help="The attributes file for info on each run [default %default]")
    )


opt  <- parse_args(
    OptionParser(#usage= "usage: %prog [options]",
        option_list=option_list)
        )

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
#

# Set column in attribute(s) file to include in DESeq objects
design <- "type"

# What groups do we want to compare?
# Which column in the attributes should be compared?
att_cols <- c("organism_part")
group1s <- c("appendix")
group2s <- c("colon")

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
                              design= ~organism_part)

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

write.table(final_results,paste(group1s[i],"_vs_",group2s[i],"results.txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
}


