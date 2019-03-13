# paste accessions into the below link
# gs://ncbi_sra_rnaseq/"<run accession>".bam

require(tidyverse)
require(magrittr)

accessions <- read_tsv("RNA-Seq-in-the-Cloud/Data/studies", col_names = FALSE)
accessions <- rename(accessions, study = X1)

accessions %<>% 
  mutate(query = paste(paste("gs://ncbi_sra_rnaseq/", accessions$study, sep = ""), ".bam", sep = ""))

library(SRAdb)
# These lines take forever. Don.t load twice.
sqlfile <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')
timeStart <- proc.time()
sqlfile <- getSRAdbFile()
proc.time() - timeStart

# connect
sra_con <- dbConnect(SQLite(),sqlfile)
dbListFields(sra_con,"study")

#convert between study and experiment ID
sraConvert(in_acc=accessions$study, out_type=c('experiment'), sra_con=sra_con ) %>% head()

#create possible cancer table
possible_cancer_samples <- getSRA(search_terms ='study_title: cancer', out_types=c('sra', 'submission', 'study', 'experiment', 'sample', 'run'), sra_con=sra_con)

possible_cancer_samples %<>% select(experiment, study_title, study_abstract, study_description, design_description, description) %>% unique()

possible_cancer_samples$normals_description <- grepl("normal", possible_cancer_samples$description, ignore.case = TRUE)
possible_cancer_samples$tumors_description <- grepl("tumor", possible_cancer_samples$description, ignore.case = TRUE)
# View(possible_cancer_samples)

hack_data <- read.table(file = "RNA-Seq-in-the-Cloud/Metadata/data/experiment_to_terms.tsv", header = FALSE, sep = "\t")

possible_cancer_samples <- unique(possible_cancer_samples$experiment)

query_disease_metadata <- subset(hack_data, possible_cancer_samples %in% hack_data$V1)

query_disease_metadata_table <- as.data.frame(table(query_disease_metadata$V2)) 

query_disease_metadata_sorted <- query_disease_metadata_table[order(query_disease_metadata_table$Freq, decreasing = TRUE),]

query_disease_metadata_sorted$percentage <- query_disease_metadata_sorted$Freq/query_disease_metadata_sorted[1,2]
query_disease_metadata_top10 <- query_disease_metadata_sorted[which(query_disease_metadata_sorted$percentage >= 0.1),]

query_disease_metadata_top10$Metadata <- query_disease_metadata_top10$Var1

ggplot(query_disease_metadata_top10, aes(Metadata, percentage)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + ggtitle("Metadata of Possible Cancer Samples")

