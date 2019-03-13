## Metadata table (available in 10 percent or more of samples)
query_disease_metadata_top10_table <- data.frame(query_disease_metadata_top10$Metadata, query_disease_metadata_top10$Freq, query_disease_metadata_top10$Proportion)
query_disease_metadata_top10_table <- setNames(query_disease_metadata_top10_table, c("Metadata", "Sample Number", "Proportion"))
query_disease_metadata_top10_table
