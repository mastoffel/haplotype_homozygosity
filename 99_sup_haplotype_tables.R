# top haplotype tables

library(snpStats)
library(tidyverse)
library(data.table)
library(here)
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

top_hap_table <- map_dfr(1:3, hap_bases) %>% 
                        mutate(snps = paste0(start_snp, "-\n", end_snp)) %>% 
                        select(-start_snp, -end_snp) %>% 
                        mutate(pos = paste0(start_pos, "-\n", end_pos)) %>% 
                        select(-start_pos, -end_pos)
library(gt)
gt(top_hap_table)

# make function to insert \n evert 100 chars
library(stringi)
ins_linebreak <- function(s) {
        paste(c(substr(s, 1, 50), "\n",
              substr(s, 51, 100), "\n", 
              substr(s, 101, 150), "\n",
              substr(s, 151, 200), "\n",
              substr(s, 201, 250), "\n", 
              substr(s, 251, 300), "\n",
              substr(s, 301, 350), "\n",
              substr(s, 351, 400)), collapse = "")
}
library(scales)
hap_table <- top_hap_table %>% 
        rowwise() %>% 
        mutate(hap = ins_linebreak(hap),
               hap_bases = ins_linebreak(hap_bases)) %>% 
        mutate(p_val = scientific(p_val, digits = 3))


hap_table1 <- select(hap_table, -hap)
hap_table2 <- select(hap_table, -hap_bases)
gt(hap_table1) %>% 
        cols_label(
                #hap = "Haplotype binary (0 = Reference / 1 = Alternative Allele)",
                hap_bases = "Haplotype",
                p_val = "p value (chi-square test)",
                obs = "# observed\nhomozyous\noffspring",
                exp = "# expected\nhomozygous\noffspring",
                num_carriers = "# carrier x carrier matings",
                chr = "chromosome",
                snps = "first SNP / last SNP",
                pos = "start position / end position (bp)"
                # start_pos = "start (bp)",
                # end_pos = "end (bp)",
                # start_snp = "start SNP",
                # end_snp = "end SNP"
                
        ) %>% 
        cols_width(
                p_val ~ px(130),
                num_carriers ~ px(100),
                pos ~ px(150),
                chr ~ px(50)
        ) %>% 
        cols_align(
                align = c("left"),
                columns = everything()
        ) %>% 
        gtsave(
                "tables/haplotypes.png",
                vwidth = 1400, vheight = 1000
        )
