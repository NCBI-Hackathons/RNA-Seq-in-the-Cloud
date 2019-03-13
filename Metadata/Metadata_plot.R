hack_data <- read.table(file = "experiment_to_terms.tsv.txt", header = FALSE, sep = "\t")
library(rjson)
query_disease_samples <- fromJSON(file = "/Users/khunzawlatt/Desktop/NCBI_hackathon/hackathon_clone/RNA-Seq-in-the-Cloud/Metadata/data/experiment_to_terms_term-breast-cancer.json")


query_disease_metadata <- subset(hack_data, V1 %in% query_disease_samples)

query_disease_metadata_table <- as.data.frame(table(query_disease_metadata$V2)) 

query_disease_metadata_sorted <- query_disease_metadata_table[order(query_disease_metadata_table$Freq, decreasing = TRUE),]
query_disease_metadata_sorted$Metadata <- query_disease_metadata_sorted$Var1

query_disease_metadata_sorted$Proportion <- query_disease_metadata_sorted$Freq/query_disease_metadata_sorted[1,2]
query_disease_metadata_top10 <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Proportion >= 0.1),]



library(ggplot2)

## Bar graph for proportion of top metadata
bp <- ggplot(query_disease_metadata_top10, aes(x = reorder(Metadata, -Proportion), y = Proportion)) +
  xlab("Metadata") + 
  ggtitle("Metadata available in 10 percent or more of samples") +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

bp
