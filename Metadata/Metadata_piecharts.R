library(ggplot2)

par(mfrow=c(2,2))

## Pie charts for top metadata of interest
query_cell_line <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "cell line"),]
query_cell_line_pie <- c(query_cell_line$Proportion, 1 - query_cell_line$Proportion)
pie_cell_line_labels <- c("cell line", "no information")
pie_cell_line <- pie(query_cell_line_pie, labels = pie_cell_line_labels, main="Cell Line", col = rainbow(length(query_cell_line_pie)))

query_treatment <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "treatment"),]
query_treatment_pie <- c(query_treatment$Proportion, 1 - query_treatment$Proportion)
pie_treatment_labels <- c("treatment", "no information")
pie_treatment <- pie(query_treatment_pie, labels = pie_treatment_labels, main="Treatment", col = rainbow(length(query_treatment_pie)))

query_adult <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "adult"),]
query_adult_pie <- c(query_adult$Proportion, 1 - query_adult$Proportion)
pie_adult_labels <- c("adult", "no information")
pie_adult <- pie(query_adult_pie, labels = pie_adult_labels, main="Adult", col = rainbow(length(query_adult_pie)))

query_female <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "female organism"),]
query_male <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "male organism"),]
query_male_female_pie <- c(query_male$Proportion, query_female$Proportion, 1 - (query_female$Proportion + query_male$Proportion))
pie_male_female_labels <- c("male", "female", "no information")
pie_male_female <- pie(query_male_female_pie, labels = pie_male_female_labels, main="Sex", col = rainbow(length(query_male_female_pie)))
