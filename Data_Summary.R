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
# sqlfile <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')
# timeStart <- proc.time()
# sqlfile <- getSRAdbFile()
# proc.time() - timeStart

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
View(possible_cancer_samples)
