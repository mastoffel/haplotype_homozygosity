# repeat genedropping with loss of haps due to recombination

library(tidyverse)
cm_map <- read_table("../sheep_roh/data/sheep_geno_imputed_oar31_17052020_cM_full.txt",
                     col_names = FALSE) %>% 
                     setNames(c("chr", "snp", "cM", "bp", "cMPosition", "org_map", "index"))

sel18_cm_length <- cm_map %>% 
                        filter(chr == 18 & (bp == 4691643 | bp == 7448353)) %>% 
                        summarise(diff(cMPosition)) %>% 
                        .[[1]]
                        
sel07_cm_length <- cm_map %>% 
        filter(chr == 7 & bp == 71164579 | bp == 73326795) %>% 
        summarise(diff(cMPosition)) %>% 
        .[[1]]

sel05_cm_length <- cm_map %>% 
        filter(chr == 5 & bp == 37164925| bp == 39808507) %>% 
        summarise(diff(cMPosition)) %>% 
        .[[1]]

