library(ggplot2)

piecharts <- par(mfrow=c(2,2))

## Pie charts for top metadata of interest
#query_cell_line <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "cell line"),]

query_cell_line <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "cell line"),]
if (nrow(query_cell_line) == 0) {
  query_cell_line_pie <- c(0.0, 1.0)
} else {
  query_cell_line_pie <- c(query_cell_line$Proportion, 1 - query_cell_line$Proportion)
}
pie_cell_line_labels <- c("cell line", "no information")
pie_cell_line <- pie(query_cell_line_pie, labels = pie_cell_line_labels, main="Cell Line", col = rainbow(length(query_cell_line_pie)))


query_treatment <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "treatment"),]
if (nrow(query_treatment) == 0) {
  query_treatment_pie <- c(0.0, 1.0)
} else {
query_treatment_pie <- c(query_treatment$Proportion, 1 - query_treatment$Proportion)
}
pie_treatment_labels <- c("treatment", "no information")
pie_treatment <- pie(query_treatment_pie, labels = pie_treatment_labels, main="Treatment", col = rainbow(length(query_treatment_pie)))

query_adult <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "adult"),]
if (nrow(query_adult) == 0) {
  query_adult_pie <- c(0.0, 1.0)
} else {
  query_adult_pie <- c(query_adult$Proportion, 1 - query_adult$Proportion) 
}
pie_adult_labels <- c("adult", "no information")
pie_adult <- pie(query_adult_pie, labels = pie_adult_labels, main="Adult", col = rainbow(length(query_adult_pie)))

query_female <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "female organism"),]
query_male <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$Metadata == "male organism"),]

if (nrow(query_male) == 0) {
  male_prop <- 0.0
} else {
  male_prop <- query_male$Proportion
}


if (nrow(query_female) == 0) {
  female_prop <- 0.0
} else {
  female_prop <- query_female$Proportion
}

query_male_female_pie <- c(male_prop, female_prop, 1 - (male_prop + female_prop))
pie_male_female_labels <- c("male", "female", "no information")
pie_male_female <- pie(query_male_female_pie, labels = pie_male_female_labels, main="Sex", col = rainbow(length(query_male_female_pie)))

piecharts
