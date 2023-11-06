library(biomaRt)
listEnsembl()

ensembl <- useEnsembl(biomart = "genes")

datasets <- listDatasets(ensembl)
head(datasets)
searchDatasets(mart = ensembl, pattern = "oaries")
ensembl <- useDataset(dataset = "oaries_gene_ensembl", mart = ensembl)

# get ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "oaries_gene_ensembl")

#
filters = listFilters(ensembl)
filters[1:5,]

attributes = listAttributes(ensembl)
attributes[1:5,]

# hap chr 9 38292279
start <- 38292279 - 250000
end <- 38292279 + 500000
# hap 18 2915687
start <- 2915687 - 500000
end <- 2915687 + 500000
# hap 5 36669706
start <- 36669706 - 500000
end <- 36669706 + 500000
# hap 7 70613419
start <- 70613419 - 250000
end <- 70613419 + 250000

res <- getBM(attributes = c("external_gene_name","ensembl_gene_id","description", "chromosome_name", "start_position", "end_position"),
      filters = c("chromosome_name", "start", "end"),
      values = list(7, start, end),
      mart = ensembl)

res %>% 
        as_tibble() %>% 
        mutate(start_position = start_position/1e6,
               end_position = end_position/1e6,
               hit = 38292279/1e6)

res <- getBM(attributes = c("external_gene_name", "description", "chromosome_name", "start_position", "end_position"),
             filters = c("chromosome_name"),
             values = list(18),
             mart = ensembl)
sum(res$external_gene_name == "U6")
