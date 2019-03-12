# Biomart
library(biomaRt)
mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# listAttributes(mart) %>% View()


#Get HGNC names from transcript_ID's and 
getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "description"),
      filters    = "hgnc_symbol",
      values     = "A1BG-AS1", 
      mart       = mart) %>% View()

# As opposed to "A1BG-AS1", insert a vector of transcript IDs. Also, remove the "description" argument rom the above to remove the functional annotation.