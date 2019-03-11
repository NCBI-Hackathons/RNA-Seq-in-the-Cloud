library(DESeq2)

# Set paths for two dataframes created by Tomer
# 1. read counts for each comparable run
counts_path <- ""
#
# 2. attributes for each run
atts_path <- ""

# Set column in attribute(s) file to include in DESeq objects
design <- 

# What groups do we want to compare?
# i.e. which column in the attributes should be compared?
att_cols <- c("")
group1s <- c("")
group2s <- c("")


counts <- read.table(counts_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
atts <- read.table(atts_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# create DESeq Object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = atts,
                              design= design)

dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 2 ) >= 2
dds <- dds[idx,]
dds <- DESeq(dds)
rld <- vst(dds)

for (i in c(1:length(att_col)))
{
att_cols[i]
group1s[i]
group2s[i]
  
# contrast groups
results(dds, contrast=c(att_cols[i],
                        group1s[i],
                        group2s[i]))

#rld.sub <- rld[ , rld$stage %in% c("8dpf") ]
plotPCA(rld, intgroup=c(att_cols[i]))+theme_bw()

pcaData <- plotPCA(rld, intgroup=c(att_cols[i]), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p1 <- ggplot(pcaData, aes(PC1, PC2, color=att_cols[i])) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme_bw()+scale_shape_manual(values=c(17, 16))
#tiff("pca.tiff", width = 6, height = 6, units = 'in', res = 1000)
p1 
#dev.off()
}