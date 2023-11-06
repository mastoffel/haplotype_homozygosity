# Genedrop analysis for haplotypes
# ->Figure 2

library(genedroppeR)
library(tidyverse)
library(snpStats)
library(data.table)
library(here)
library(glue)

# desired format
data(unicorn)
unicorn %>% tibble()

# pedigree
load(here("data", "sheep_ped.RData"))
ped <- as_tibble(sheep_ped) %>% 
        setNames(c("id", "mother", "father"))

# top haplotypes (genome-wide sig. homozygote deficiency)
top_haps <- read_delim(here("output", "top_haps_400.txt"))

# individual haplotypes and fitness
haps_fit <- read_delim(here("output", "haps400_and_fitness.txt")) %>% 
                filter(age == 0)

birth_years <- read_delim("../sheep_ID/data/1_Annual_Fitness_Measures_April_20190501.txt", 
                          delim = "\t") %>% 
        select(ID, BIRTHYEAR) %>% 
        rename(id = ID, birth_year = BIRTHYEAR) %>% 
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
                rename(cohort = birth_year) %>% 
                filter((cohort >= 1990)) %>% 
                filter(!(id == "2319" & mother == "2232" & father == "M0065"),
                       !(id == "5586" & mother == "5635" & father == "5117"),
                       !(id == "6210" & mother == "6780" & father == "M0066"),
                       !(id == "41" & father == "2332"),
                       !(id == "9745" & father == "8701")) %>% 
                mutate(gt = ifelse(gt == "AB", "BA", gt)) %>% 
                mutate(gt = as.factor(gt))
        
        # drop
        sheep_UF <- genedrop_snp(id = sheep$id,
                                 mother = sheep$mother,
                                 father = sheep$father,
                                 cohort = sheep$cohort,
                                 genotype = sheep$gt,
                                 nsim = 1000,
                                 n_founder_cohorts = 3, #8
                                 fix_founders = T,
                                 verbose = T,
                                 interval = 50)
        
        write_delim(sheep_UF, here("output", glue("genedrop3_", hap, ".txt")))
        
}

haps <- c("chr5_6293","chr7_12196","chr18_267")
walk(haps, gene_drop_hap)

