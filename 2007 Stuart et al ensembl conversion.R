library(dplyr)

ensembl_input <- read.csv("../../stuart_genes.csv") 

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

output <- data.frame()

output <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene", "description"),
  values = ensembl_input,
  mart = mart
)

write.csv(output, file = "stuart_output.csv")
listAttributes(mart)
