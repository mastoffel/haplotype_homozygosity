# top haplotype tables

library(snpStats)
library(tidyverse)
library(data.table)

top_haps <- read_delim(here("output", "top_haps_400.txt"))
top_haps

# snp map
snp_map <- fread(here("data", "plink", "sheep.bim")) %>% 
        setNames(c("chr", "snp", "cM", "pos", "a", "b")) %>% 
        #select(chr, snp, pos) %>% 
        group_by(chr) %>% 
        mutate(number = 1,
               snp_num = cumsum(number)) %>% 
        select(-number)

snp_map

top_haps

hap_bases <- function(top_hap_num){
        
        top_hap <- top_haps[top_hap_num, ]
        
        snp_map_mod <- snp_map %>% 
                filter(chr == top_hap$chr) %>% 
                filter(snp_num >= top_hap$snp_start & snp_num < top_hap$snp_start + 400) %>% 
                mutate(hap =strsplit(top_hap$hap[[1]], split = "")[[1]]) %>% 
                mutate(hap_bases = ifelse(hap == 0, a, b))
        
        top_hap %>% 
                mutate(hap_bases = paste(snp_map_mod$hap_bases, collapse = "")) %>% 
                mutate(start_snp = snp_map_mod[1, "snp"][[1]],
                       start_pos = snp_map_mod[1, "pos"][[1]],
                       end_snp = snp_map_mod[nrow(snp_map_mod), "snp"][[1]],
                       end_pos = snp_map_mod[nrow(snp_map_mod), "pos"][[1]]) %>% 
                select(hap, hap_bases, p_val, obs, exp, num_carriers, chr, start_pos, end_pos, start_snp, end_snp)
        
}

top_hap_table <- map_dfr(1:3, hap_bases)

