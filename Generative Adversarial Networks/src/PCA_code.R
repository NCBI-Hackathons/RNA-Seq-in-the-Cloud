# Combined PCA analysis
df_normalized_case <- read.csv("normalized_case.csv", header = FALSE);
df_normalized_control <- read.csv("normalized_control.csv", header = FALSE);
df_raw_case <- read.csv("raw_case.csv", header = FALSE);
df_raw_control <- read.csv("raw_control.csv", header = FALSE);

# View dimentions of the files
dim(df_normalized_case)
dim(df_normalized_control)
dim(df_raw_case)
dim(df_raw_control)

# Makes everything numeric
df_normalized_case_2 <- as.data.frame(sapply(df_normalized_case, as.numeric))
df_normalized_control_2 <- as.data.frame(sapply(df_normalized_control, as.numeric))
df_raw_case_2 <- as.data.frame(sapply(df_raw_case, as.numeric))
df_raw_control_2 <- as.data.frame(sapply(df_raw_control, as.numeric))

# Add index column
df_normalized_case_2$newcolumn <- 1
df_normalized_control_2$newcolumn <- 2
df_raw_case_2$newcolumn <- 3
df_raw_control_2$newcolumn <- 4


# Combine dataframes
combined_dataset <- rbind(df_normalized_case_2, df_normalized_control_2) # normalized data

combined_dataset_2 <- rbind(df_raw_case_2, df_raw_control_2) # raw data

# Actual PCA
library(ggfortify)
autoplot(prcomp(combined_dataset), colour = 'newcolumn')
autoplot(prcomp(combined_dataset_2), colour = 'newcolumn')

# Get PCA information
library("factoextra")

# Extract the results for variables and individuals
result_normalized <- get_pca(combined_dataset)
result_raw <- get_pca(combined_dataset_2)

write.csv(result_normalized, "result_normalized.csv")
write.csv(result_raw, "result_raw.csv")







