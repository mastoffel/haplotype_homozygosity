library(tidyverse)
library(data.table)
library(here)
library(snpStats)

# snp map
snp_map <- fread(here("data", "plink", "sheep.bim")) %>% 
        setNames(c("chr", "snp", "cM", "pos", "a", "b")) %>% 
        select(chr, snp, pos) %>% 
        group_by(chr) %>% 
        mutate(number = 1,
               snp_num = cumsum(number)) %>% 
        select(-number)


# find out which snp differs
top_haps <-  read_delim(here("output", "top_haps_500.txt"), col_types = (c("cdccccddddcc")))
start <- top_haps %>% filter(chr == 9) %>% .$snp_start
hap1 <- str_split(top_haps$haps[1], "_")[[1]][1]
hap2 <- str_split(top_haps$haps[1], "_")[[1]][2]

snp_num <- which(!(str_split(hap1, "")[[1]] == str_split(hap2, "")[[1]]))
snp_ind <- start + snp_num - 1

snp <- snp_map %>% 
        filter(chr == 9) %>% 
        filter(snp_num == snp_ind)


# no individual has less than 95 %
sheep_plink_name <- "../sheep/data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2_snpfile_for_merge"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

snp_geno <- as_tibble(as(full_sample$genotypes[, snp$snp], Class = "numeric"),
                                  rownames = "id")

sheep_plink_name <- "../sheep/data/SNP_chip/oar31_mapping/Plates_1to87_QC3_snpfile_for_merge"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample_50k <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snp_geno_50k

snp_geno_50k <- as_tibble(as(full_sample_50k$genotypes[, snp$snp], Class = "numeric"),
                      rownames = "id") %>% 
                rename(s23535.1_50K=s23535.1)

compare <- snp_geno %>% 
        left_join(snp_geno_50k) %>% 
        mutate(same = s23535.1 == s23535.1_50K) 
sum(compare$same)
