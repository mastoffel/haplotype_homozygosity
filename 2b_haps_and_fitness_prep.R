library(dplyr)
library(readr)
library(here)
library(glue)
library(tidyr)
library(purrr)
library(tidyverse)
library(stringr)
library(data.table)
library(broom.mixed)
library(wesanderson)
source("theme_simple.R")


# fitness
load("data/fitness_roh.RData") 
age <- fitness_data %>% 
        select(id, sheep_year, age) %>% 
        rename(mum_age = age)

fitness <- fitness_data %>% 
                select(id, survival, sheep_year, age, birth_year, sex, mum_id,
                       twin, offspring_born, offspring_survived, weight,
                       hindleg, horn, father, froh_all, birthwt) %>% 
                group_by(id) %>% 
                mutate(lifespan = max(age)) %>% 
                # mother age
                left_join(age, by = c("mum_id" = "id", "sheep_year"))

# snp map
snp_map <- fread(here("data", "plink", "sheep.bim")) %>% 
  setNames(c("chr", "snp", "cM", "pos", "a", "b")) %>% 
  select(chr, snp, pos) %>% 
  group_by(chr) %>% 
  mutate(number = 1,
         snp_num = cumsum(number)) %>% 
  select(-number)

# target haplotypes:
# read results from haplotype homozygosity scan
all_files <- list.files(here("output", "hap_results_imputed", "hap_len_400"), full.names = TRUE)
results <- map(all_files, read_delim, delim = "\t") %>%  #read_delim, delim = "\t"
        bind_rows()

# note: sometimes there are two significant haplotypes starting with the same 
# snp -> need to be merged eventually
top_haps <- results %>% 
        # get haplotypes which are significant because there are fewer homozyous
        # haplotypes as expected
        filter(obs < exp) %>% 
        filter(p_val < (0.05/39184)) %>% # 417373 # 39184
        arrange(chr, snp_start) %>% 
        # workflow to get only one haplotype per genome-region
        group_by(chr) %>% 
        # work with lag, the default gets around the a problem
        # when the haplotype starts at the very beginning of a chromosome
        # (snp numbers below 10)
        mutate(lag1 = lag(snp_start),
               diff = snp_start - lag1,
               diff = ifelse(is.na(diff), snp_start, diff)) %>% 
        mutate(region_start = ifelse(diff < 100, NA, snp_start)) %>% 
        fill(region_start) %>% 
        mutate(region = paste0("chr", chr, "_", region_start)) %>% 
        group_by(region, snp_start) %>% 
        # if multiple haplotypes with same starting snp are significant, 
        # add them to "haps" column
        mutate(haps = paste0(hap, collapse = "_"),
               obs = paste0(obs, collapse = ","),
               exp = paste0(exp, collapse = ","),
               num_carriers = paste0(num_carriers, collapse = ",")) %>% 
        # now haplotypes with lowest pvalues per regions
        group_by(region) %>% 
        arrange(p_val) %>% 
        filter(row_number()==1) 
      
top_haps
write_delim(top_haps, here("output", "top_haps_400.txt"))
top_haps <- read_delim(here("output", "top_haps_400.txt"))
# on chromosome 9, it's two significant haplotypes at the same location
# they only differ by one mutation
# hap_seq <- unlist(str_split(top_haps[1, ]$haps, "_"))
# diff_haps <- strsplit(hap_seq[1], "")[[1]] == strsplit(hap_seq[2], "")[[1]]
# which(!diff_haps)
# strsplit(hap_seq[2], "")[[1]][!diff_haps]
# snp_map %>% 
#   filter(chr == 14) %>% 
#   filter(snp_num == (6599+ which(!diff_haps) - 1))

# take top haplotype and get genotypes for hom/het/alt_hom 
hap_to_geno <- function(i, top_haps) {
        
        top_hap <- top_haps[i, ]
        
        chr <- top_hap$chr
        
        # get haplotypes, remove everything but genotypes to make it a matrix
        haps_raw <- fread(here("output", "phased_matrix", glue("sheep_phased_chr_{chr}.hap.gz"))) %>% 
                select(-V1, -V2, -V3, -V4, -V5) %>% 
                as.matrix()
        
        # get individual names, should be 5952
        # here is a bit clunky here so stay with fixed path
        ind_names <- system(glue("zgrep '^#CHROM*' output/phased/sheep_phased_chr_{chr}.vcf.gz"),
                            intern = TRUE) %>% 
                str_split("\t") %>%
                unlist() %>% 
                .[-(1:9)] %>% 
                str_split("_") %>% 
                map_chr(2)
        
        # double each name and add _a _b for haplotype 1 / haplotype 2
        ind_names_hap <- rep(ind_names, each = 2) %>% 
                paste0(., c("_a", "_b"))
        # add to matrix
        colnames(haps_raw) <- ind_names_hap
        
        snp_start <- top_hap$snp_start
        hap_length <- str_count(top_hap$hap)
        
        # combine genotypes to haplotypes
        haps <- apply(haps_raw[snp_start:(snp_start+hap_length-1), ],
                      paste, collapse = "", MARGIN = 2) %>% 
                enframe(name = "id_hap", value = "hap") %>% 
                mutate(id = str_remove(id_hap, "_[a-z]")) %>%
                select(-id_hap) %>% 
                mutate(hap_pos = rep(c("hap_a", "hap_b"), nrow(.)/2), .before = 1) %>% 
                mutate(id = as.numeric(id)) %>% 
                pivot_wider(names_from = hap_pos, values_from = hap)
        
        # in one position, two haplotypes were significant
        hap_seq <- unlist(str_split(top_hap$haps, "_"))
        # genotype for haplotype
        # in case of multiple significant haplotypes, this assumes that all
        # significant backgrounds contain the target mutation, and are 
        # therefore clustered as one haplotype
       
        hap_gts <- haps  %>% 
          rowwise() %>% 
          mutate(gt = case_when(
            (hap_a %in% hap_seq) & (hap_b %in% hap_seq) ~ 2,
            (hap_a %in% hap_seq) & !(hap_b %in% hap_seq) ~ 1,
            !(hap_a %in% hap_seq) & (hap_b %in% hap_seq) ~ 1,
            !(hap_a %in% hap_seq) & !(hap_b %in% hap_seq) ~ 0,
            TRUE ~ NA_real_
          )) %>% 
          ungroup() %>% 
          mutate(snp_start = snp_start, chr = chr)
        
        # this is specific for the case of length(hap_seq) == 2, i.e. 2
        # significant haplotypes at the same location
        if (length(hap_seq) > 1){
          hap_gts$gt_hap2 <- rowSums(hap_seq[2] == as.matrix(hap_gts[c('hap_a', 'hap_b')]))
        } else {
          hap_gts$gt_hap2 <- NA
        }
        
        hap_gts
        
   
}

# get haplotype genotypes for top haplotypes
haps <- map(1:nrow(top_haps), hap_to_geno, top_haps)
#saveRDS(haps, file="output/haps_600.rds")

# data.frame with all haplotypes and regions
haps_all <- haps %>% 
        setNames(top_haps$region) %>% 
        bind_rows(.id = "region") %>% 
        mutate(id = as.character(id)) %>% 
        left_join(fitness) #%>% 
        #mutate(gt = as.factor(gt))

haps_all %>% 
        write_delim(here("output", "haps400_and_fitness.txt"), " ")
# top haplotype genotypes per individual
# save
haps_ind <- haps %>% 
            setNames(top_haps$region) %>% 
            bind_rows(.id = "region") %>% 
            write_delim(here("output", "sheep_top_haps_400.txt"), " ")



# haplotype frequency over time
#write_delim(haps_all, here("output", "haps_all_500.txt"), " ")

haps_all %>% 
  filter(as.numeric(as.character(birth_year)) > 1990) %>% 
  group_by(birth_year, region) %>% 
  summarise(freq = sum(gt) / (n()*2)) %>% 
  ungroup() %>% 
  ggplot(aes(birth_year, freq, group = 1)) +
    geom_line() +
    facet_wrap(~region, scales = "free_y")

# check the two haps on 9
haps_all %>% 
  filter(as.numeric(as.character(birth_year)) > 1990) %>% 
  filter(chr == 9 & gt_hap2 == 1) %>% 
  group_by(birth_year) %>% 
  summarise(freq = sum(gt) / (n()*2)) %>% 
  ggplot(aes(birth_year, freq, group = 1)) +
  geom_line()
