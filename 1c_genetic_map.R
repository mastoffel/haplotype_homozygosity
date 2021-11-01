library(tidyverse)
library(glue)
# make genetic map per chromosome for shapeit4
# structure is:

# pos   chr     cM
# 10894949      21      0.052423

# interpolated map from "Mutation load etc. paper"
# files: sheep.bim (imputed) \ sheep_50K.bim
plink_map <- read_delim("data/plink/sheep_50K.bim", "\t", col_names = FALSE) %>%  
                select(X2, X4) %>% 
                setNames(c("snp", "pos"))
genetic_map <- read_delim("../sheep_roh/data/sheep_geno_imputed_oar31_17052020_cM_full.txt", # for imputed data
                          delim = "\t", col_names = FALSE) %>% 
                                select(X2, X4, X1, X5) %>% 
                                setNames(c("snp", "pos", "chr", "cM")) %>% 
                          # join snps present in plink files
                          right_join(plink_map)
# 50K
genetic_map_50K <- read_delim("../sheep_roh/data/7_20200504_Full_Linkage_Map.txt",
                              delim = "\t") %>% 
                        select(SNP.Name, cMPosition, Chr) %>% 
                        setNames(c("snp", "cM", "chr")) %>% 
                        right_join(plink_map, by = "snp")

#write_delim(genetic_map, "data/genetic_map_all_chr.txt", "\t")

# split by chromosome and write as gz

write_by_chr <- function(chr_num, genetic_map) {
        genetic_map %>% 
                select(-snp) %>% 
                filter(chr == chr_num) %>% 
                write_delim(glue("data/genetic_map_50K/chr_{chr_num}.gmap.gz"))
}

walk(1:26, write_by_chr, genetic_map)




