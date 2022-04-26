# genedrop haplotypes
library(genedroppeR)
library(tidyverse)
library(snpStats)
library(data.table)
library(here)

# desired format
data(unicorn)
unicorn %>% tibble()

# pedigree
load(here("data", "sheep_ped.RData"))
ped <- as_tibble(sheep_ped) %>% 
        setNames(c("id", "mother", "father"))

# top_haplotypes
top_haps <- read_delim(here("output", "top_haps_400.txt"))
hap1 <- top_haps$hap[1]

# haplotype "chr18_267"
haps_fit <- read_delim(here("output", "haps400_and_fitness.txt")) %>% 
                filter(age == 0) %>% 
                filter(region == "chr18_267") %>%  # chr7_12196 chr18_267
                select(id, hap_a, hap_b, birth_year, sex) %>% 
                mutate(hap_a = ifelse(hap_a == hap1, "A", "B"),
                       hap_b = ifelse(hap_b == hap1, "A", "B")) %>% 
                mutate(gt = paste0(hap_a, hap_b)) %>% 
                mutate(id = as.character(id))

sheep <- ped %>% 
        left_join(haps_fit) %>% 
        select(id, mother, father, birth_year, gt) %>% 
        rename(cohort = birth_year) %>% 
        filter(cohort >= 1990) %>% 
        filter(!(id == "2319" & mother == "2232" & father == "M0065"),
               !(id == "5586" & mother == "5635" & father == "5117")) %>% 
        mutate(gt = ifelse(gt == "AB", "BA", gt)) %>% 
        mutate(gt = as.factor(gt))

sheep_summ <- summary_cohort(id = sheep$id,
               mother = sheep$mother,
               father = sheep$father,
               cohort = sheep$cohort,
               genotype = sheep$gt)

ggplot(sheep_summ, aes(cohort, A)) +
        geom_line() +
        stat_smooth(method = "lm") +
        ggtitle("Temporal dynamics of A allele")

ggplot(sheep_summ, aes(cohort, PropFounders)) +
        geom_line() +
        ggtitle("Proportion of Founder IDs per cohort")

ggplot(sheep_summ, aes(cohort, PropGenotyped)) +
        geom_line()  +
        ggtitle("Proportion of Genotyped IDs per cohort")

sheep_UF <- genedrop_snp(id = sheep$id,
                         mother = sheep$mother,
                         father = sheep$father,
                         cohort = sheep$cohort,
                         genotype = sheep$gt,
                         nsim = 500,
                         n_founder_cohorts = 8, #8
                         fix_founders = T,
                         verbose = T,
                         interval = 50)
tibble(sheep_UF)
sheep_UF_summ <- summary_genedrop(sheep_UF)
str(sheep_UF_summ)
plot_genedrop_results(sheep_UF_summ)

plot_genedrop_lm_slopes(sheep_UF_summ,
                        n_founder_cohorts = 10,
                        remove_founders = F)
plot_genedrop_cumulative_change(sheep_UF_summ)

write_delim(sheep_UF, here("output", "genedrop_hap_chr18_267.txt"))


# haplotype "chr18_267"
hap1 <- top_haps$hap[2]
haps_fit <- read_delim(here("output", "haps400_and_fitness.txt")) %>% 
        filter(age == 0) %>% 
        filter(region == "chr7_12196") %>% 
        select(id, hap_a, hap_b, birth_year, sex) %>% 
        mutate(hap_a = ifelse(hap_a == hap1, "A", "B"),
               hap_b = ifelse(hap_b == hap1, "A", "B")) %>% 
        mutate(gt = paste0(hap_a, hap_b)) %>% 
        mutate(id = as.character(id))

sheep <- ped %>% 
        left_join(haps_fit) %>% 
        select(id, mother, father, birth_year, gt) %>% 
        rename(cohort = birth_year) %>% 
        filter(cohort >= 1990) %>% 
        filter(!(id == "2319" & mother == "2232" & father == "M0065"),
               !(id == "5586" & mother == "5635" & father == "5117")) %>% 
        mutate(gt = ifelse(gt == "AB", "BA", gt)) %>% 
        mutate(gt = as.factor(gt))

sheep_summ <- summary_cohort(id = sheep$id,
                             mother = sheep$mother,
                             father = sheep$father,
                             cohort = sheep$cohort,
                             genotype = sheep$gt)

ggplot(sheep_summ, aes(cohort, A)) +
        geom_line() +
        stat_smooth(method = "lm") +
        ggtitle("Temporal dynamics of A allele")

ggplot(sheep_summ, aes(cohort, PropFounders)) +
        geom_line() +
        ggtitle("Proportion of Founder IDs per cohort")

ggplot(sheep_summ, aes(cohort, PropGenotyped)) +
        geom_line()  +
        ggtitle("Proportion of Genotyped IDs per cohort")

sheep_UF <- genedrop_snp(id = sheep$id,
                         mother = sheep$mother,
                         father = sheep$father,
                         cohort = sheep$cohort,
                         genotype = sheep$gt,
                         nsim = 500,
                         n_founder_cohorts = 10,
                         fix_founders = T,
                         verbose = T,
                         interval = 100)
tibble(sheep_UF)
sheep_UF_summ <- summary_genedrop(sheep_UF)
str(sheep_UF_summ)
plot_genedrop_results(sheep_UF_summ)

plot_genedrop_lm_slopes(sheep_UF_summ,
                        n_founder_cohorts = 10,
                        remove_founders = F)

write_delim(sheep_UF, here("output", "genedrop_hap_chr7_12196.txt"))



# haplotype "chr5_6293"
hap1 <- top_haps$hap[3]
haps_fit <- read_delim(here("output", "haps400_and_fitness.txt")) %>% 
        filter(age == 0) %>% 
        filter(region == "chr5_6293") %>% 
        select(id, hap_a, hap_b, birth_year, sex) %>% 
        mutate(hap_a = ifelse(hap_a == hap1, "A", "B"),
               hap_b = ifelse(hap_b == hap1, "A", "B")) %>% 
        mutate(gt = paste0(hap_a, hap_b)) %>% 
        mutate(id = as.character(id))

sheep <- ped %>% 
        left_join(haps_fit) %>% 
        select(id, mother, father, birth_year, gt) %>% 
        rename(cohort = birth_year) %>% 
        filter(cohort >= 1990) %>% 
        filter(!(id == "2319" & mother == "2232" & father == "M0065"),
               !(id == "5586" & mother == "5635" & father == "5117")) %>% 
        mutate(gt = ifelse(gt == "AB", "BA", gt)) %>% 
        mutate(gt = as.factor(gt))

sheep_summ <- summary_cohort(id = sheep$id,
                             mother = sheep$mother,
                             father = sheep$father,
                             cohort = sheep$cohort,
                             genotype = sheep$gt)

ggplot(sheep_summ, aes(cohort, A)) +
        geom_line() +
        stat_smooth(method = "lm") +
        ggtitle("Temporal dynamics of A allele")

ggplot(sheep_summ, aes(cohort, PropFounders)) +
        geom_line() +
        ggtitle("Proportion of Founder IDs per cohort")

ggplot(sheep_summ, aes(cohort, PropGenotyped)) +
        geom_line()  +
        ggtitle("Proportion of Genotyped IDs per cohort")

sheep_UF <- genedrop_snp(id = sheep$id,
                         mother = sheep$mother,
                         father = sheep$father,
                         cohort = sheep$cohort,
                         genotype = sheep$gt,
                         nsim = 500,
                         n_founder_cohorts = 10,
                         fix_founders = F,
                         verbose = T,
                         interval = 100)
tibble(sheep_UF)
sheep_UF_summ <- summary_genedrop(sheep_UF)
str(sheep_UF_summ)
plot_genedrop_results(sheep_UF_summ)

write_delim(sheep_UF, here("output", "genedrop_hap_chr5_6293.txt"))


