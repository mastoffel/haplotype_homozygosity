library(genedroppeR)
library(tidyverse)
library(snpStats)
library(data.table)
library(here)
library(glue)
source("theme_simple.R")

# data
haps_all <- read_delim(here("output", "haps400_and_fitness.txt"))

# top haplotypes
hap_names <- c("chr18_267", "chr7_12196", "chr5_6293")

# load genedropping simulations for all three haplotypes
load_gd <- function(hap_name) {
        read_delim(here("output", paste0("genedrop_", hap_name, ".txt")))
}
genedrops <- map(hap_names, load_gd)

# dfs for observed and summaries
gd_observed <- map(genedrops, function(x) summary_genedrop(x)$observed_frequencies)
names(gd_simulated) <- hap_names
gd_simulated <- map(genedrops, function(x) summary_genedrop(x)$simulated_frequencies)
names(gd_simulated) <- hap_names

# coerce to single df
gds_df <- bind_rows(gd_simulated, .id = "name") %>% 
                as_tibble() %>% 
                rename(region = name, birth_year = Cohort, freq = Count) %>% 
                mutate(region = factor(region, levels = c("chr5_6293", "chr7_12196",
                                                        "chr18_267")))

# colors and labels for the plot
cols <- c("#DE9151", "#F34213", "#2E2E3A")
hap_labels <- c(
        "chr5_6293" = "SEL05",
        "chr7_12196" = "SEL07",
        "chr18_267" = "SEL18"
)

p_freq <- haps_all %>% 
        mutate(region = factor(region),
               region = fct_reorder(region, chr)) %>% 
        filter(as.numeric(as.character(birth_year)) >= 1990) %>% 
        group_by(birth_year, region) %>% 
        summarise(freq = sum(gt) / (n()*2)) %>% 
        ungroup() %>% 
        #mutate(highlight = factor(ifelse(region == "chr7_12119", 1, 0))) %>% 
        ggplot(aes(birth_year, freq, group = 1, color = region)) +
        geom_line(data=gds_df, alpha = 0.7, aes(y = p, group = Simulation),
                  size = 0.02, color = "#4C566A")+
        geom_line(size = 1) +
        scale_color_manual(values = cols) +
        #geom_smooth(method = "lm", se = FALSE) +
        facet_grid(~region, labeller=labeller(region = hap_labels)) + 
        theme_simple(grid_lines = TRUE, axis_lines = TRUE) +
        scale_y_continuous(limits = c(0, 0.5)) +
        # scale_color_manual(values = c( "#5E81AC")) + #5E81AC "#94350b"
        ylab("Haplotype\nfrequency") +
        xlab("Birth year") +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.spacing = unit(2.5, "lines"),
              strip.text.x = element_blank())
# panel.grid.major.y = element_line(size = 0.1, color = "#4C566A")) 
p_freq

# 

gmap_dfr(hap_names_full, load_gd)

genedrops_sum <- summary_genedrop(genedrops)

plot_genedrop_lm_slopes(unicorn.UF.summ,
                        n_founder_cohorts = 5,
                        remove_founders = T)

















hap1 <- top_haps$hap[2]

# haplotype "chr18_267"
haps_fit <- read_delim(here("output", "haps400_and_fitness.txt")) %>% 
        filter(age == 0) %>% 
        filter(region == "chr7_12196") %>%  # chr7_12196 chr18_267 chr5_6293
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
                         nsim = 300,
                         n_founder_cohorts = 5, #8
                         fix_founders = T,
                         verbose = T,
                         interval = 50)

sheep_UF_summ <- summary_genedrop(sheep_UF)
plot_genedrop_results(sheep_UF_summ)

plot_genedrop_lm_slopes(sheep_UF_summ,
                        n_founder_cohorts = 5,
                        remove_founders = T)
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
                         nsim = 1000,
                         n_founder_cohorts = 7,
                         fix_founders = T,
                         verbose = T,
                         interval = 100)
tibble(sheep_UF)
sheep_UF_summ <- summary_genedrop(sheep_UF)
str(sheep_UF_summ)
plot_genedrop_results(sheep_UF_summ)

plot_genedrop_lm_slopes(sheep_UF_summ,
                        n_founder_cohorts = 7,
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
                         nsim = 1000,
                         n_founder_cohorts = 7,
                         fix_founders = T,
                         verbose = T,
                         interval = 100)
tibble(sheep_UF)
sheep_UF_summ <- summary_genedrop(sheep_UF)
str(sheep_UF_summ)
plot_genedrop_results(sheep_UF_summ)

write_delim(sheep_UF, here("output", "genedrop_hap_chr5_6293.txt"))


