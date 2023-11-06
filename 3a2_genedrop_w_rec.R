# Genedrop analysis for haplotypes with recombination
# -> Supplementary Figure 2

library(genedroppeR)
source("genedrop2.R")
library(dplyr)
library(furrr)
library(purrr)
library(readr)
library(snpStats)
library(data.table)
library(glue)
library(plyr)
library(here)

set.seed(32423)
# set.seed(43432)
# desired format
data(unicorn)
unicorn %>% tibble()

# pedigree
load(here::here("data", "sheep_ped.RData"))
ped <- as_tibble(sheep_ped) %>% 
        setNames(c("id", "mother", "father"))

# top haplotypes (genome-wide sig. homozygote deficiency)
top_haps <- read_delim(here::here("output", "top_haps_400.txt"))

# individual haplotypes and fitness
haps_fit <- read_delim(here::here("output", "haps400_and_fitness.txt")) %>% 
        filter(age == 0)

birth_years <- read_delim("../sheep_ID/data/1_Annual_Fitness_Measures_April_20190501.txt", 
                          delim = "\t") %>% 
        dplyr::select(ID, BIRTHYEAR) %>% 
        dplyr::rename(id = ID, birth_year = BIRTHYEAR) %>% 
        mutate(id = as.character(id)) %>% 
        group_by(id) %>% 
        sample_n(1)

# chr7_12196 chr18_267 chr5_6293
gene_drop_hap <- function(hap) {
        
        target_hap <- top_haps %>% 
                filter(region == !!hap) %>% 
                .[["hap"]]
        
        # get haplotype "genotype" for each individual
        hap_fit <- haps_fit %>% 
                filter(region == !!hap) %>%  
                select(id, hap_a, hap_b, birth_year, sex) %>% 
                mutate(hap_a = ifelse(hap_a == target_hap, "A", "B"),
                       hap_b = ifelse(hap_b == target_hap, "A", "B")) %>% 
                mutate(gt = paste0(hap_a, hap_b)) %>% 
                mutate(id = as.character(id))
        
        # combine with pedigree
        sheep <- ped %>% 
                left_join(birth_years) %>% 
                left_join(hap_fit) %>% 
                select(id, mother, father, birth_year, gt) %>% 
                dplyr::rename(cohort = birth_year) %>% 
                filter((cohort >= 1990)) %>% 
                filter(!(id == "2319" & mother == "2232" & father == "M0065"),
                       !(id == "5586" & mother == "5635" & father == "5117"),
                       !(id == "6210" & mother == "6780" & father == "M0066"),
                       !(id == "41" & father == "2332"),
                       !(id == "9745" & father == "8701")) %>% 
                mutate(gt = ifelse(gt == "AB", "BA", gt)) %>% 
                mutate(gt = as.factor(gt))
        
        # drop
        p_rec <- case_when(
                hap == "chr5_6293" ~ 0.02130331,
                hap == "chr7_12196" ~ 0.01893693,
                hap == "chr18_267" ~ 0.04415,
        )
        sheep_UF <- genedrop2(id = sheep$id,
                                 mother = sheep$mother,
                                 father = sheep$father,
                                 cohort = sheep$cohort,
                                 genotype = sheep$gt,
                                 nsim = 1000,
                                 n_founder_cohorts = 3, #8
                                 fix_founders = T,
                                 verbose = T,
                                 interval = 50,
                                 p_rec = p_rec)
        
        write_delim(sheep_UF, here("output", glue("genedrop_rec6_", hap, ".txt")))
        
}

haps <- c("chr5_6293","chr7_12196","chr18_267")

plan(multisession, workers = 3)
future_walk(haps, gene_drop_hap)

