# genedrop haplotypes
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
                select(ID, BIRTHYEAR, MOTHER, FATHER) %>% 
                rename(id = ID, birth_year = BIRTHYEAR, mother = MOTHER,
                       father = FATHER) %>% 
                mutate(id = as.character(id),
                       mother = as.character(mother),
                       father = as.character(father)) %>% 
                group_by(id) %>% 
                sample_n(1) %>% 
                select(id,birth_year) #%>% 
                #rename(birth_year_full = birth_year)

target_hap <- top_haps$hap[[1]]
hap <- top_haps$region[[1]]

# get haplotype "genotype" for each individual
hap_fit <- haps_fit %>% 
        filter(region == !!hap) %>%  
        select(id, hap_a, hap_b, sex) %>% 
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
        filter((cohort >= 1989) & (cohort < 2018)) %>% 
        filter(!(id == "2319" & mother == "2232" & father == "M0065"),
               !(id == "5586" & mother == "5635" & father == "5117"),
               !(id == "6210" & mother == "6780" & father == "M0066"),
               !(id == "41" & father == "2332"),
               !(id == "9745" & father == "8701")) %>% 
        mutate(gt = ifelse(gt == "AB", "BA", gt)) %>% 
        mutate(gt = as.factor(gt))

sheep_sum <- summary_cohort(id = sheep$id,
                               mother = sheep$mother,
                               father = sheep$father,
                               cohort = sheep$cohort,
                               genotype = sheep$gt)

ggplot(sheep_sum, aes(cohort, PropFounders)) +
        geom_line() +
        ggtitle("Proportion of Founder IDs per cohort")

ggplot(sheep_sum, aes(cohort, PropGenotyped)) +
        geom_line()  +
        ggtitle("Proportion of Genotyped IDs per cohort")


sheep_UF <- genedrop_snp(id = sheep$id,
                         mother = sheep$mother,
                         father = sheep$father,
                         cohort = sheep$cohort,
                         genotype = sheep$gt,
                         nsim = 300,
                         n_founder_cohorts = 5, #8
                         fix_founders = T,
                         verbose = T,
                         interval = 50)

sheep_UF_sum <- summary_genedrop(sheep_UF)

plot_genedrop_lm_slopes(sheep_UF_sum,
                        n_founder_cohorts = 5,
                        remove_founders = F)

plot_genedrop_cumulative_change(sheep_UF_sum,
                                n_founder_cohorts = 5,
                                remove_founders = F)
