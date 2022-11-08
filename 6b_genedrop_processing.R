library(genedroppeR)
library(tidyverse)
library(snpStats)
library(data.table)
library(here)
library(glue)
source("theme_simple.R")
source("gd_change.R")
source("gd_slopes.R")

# data
haps_all <- read_delim(here("output", "haps400_and_fitness.txt")) %>% 
                rename(hap = region)

# top haplotypes
hap_names <- c("chr18_267", "chr7_12196", "chr5_6293")

# load genedropping simulations for all three haplotypes
load_gd <- function(hap_name) {
        read_delim(here("output", paste0("genedrop_3fc_", hap_name, ".txt")),
                   col_types = "ddddddd") %>% 
                as.data.frame()
}
genedrops <- map(hap_names, load_gd)
gd_summaries <- map(genedrops, summary_genedrop)
names(gd_summaries) <- hap_names

# dfs for observed and simulations
gd_observed <- map(gd_summaries, "observed_frequencies")
names(gd_observed) <- hap_names
gd_simulated <- map(gd_summaries, "simulated_frequencies")
names(gd_simulated) <- hap_names

# coerce to single df
gds_df <- bind_rows(gd_simulated, .id = "name") %>% 
                as_tibble() %>% 
                rename(hap = name, birth_year = Cohort, freq = Count) %>% 
                mutate(hap = factor(hap, levels = c("chr5_6293", "chr7_12196",
                                                        "chr18_267")))

# colors and labels for the plot
cols <- c("#DE9151", "#F34213", "#2E2E3A")
cols <- c('#161032', '#C73E1D', "#FFC53A")
hap_labels <- c(
        "chr5_6293" = "SEL05",
        "chr7_12196" = "SEL07",
        "chr18_267" = "SEL18"
)

p_freq <- haps_all %>% 
        mutate(hap = factor(hap),
               hap = fct_reorder(hap, chr)) %>% 
        filter(as.numeric(as.character(birth_year)) >= 1990) %>% 
        group_by(birth_year, hap) %>%  
        summarise(freq = sum(gt) / (n()*2)) %>% 
        ungroup() %>% 
        #mutate(highlight = factor(ifelse(region == "chr7_12119", 1, 0))) %>% 
        ggplot(aes(birth_year, freq, group = 1, color = hap)) +
        geom_line(data=gds_df, 
                  alpha = 1, aes(y = p, group = Simulation),
                  size = 0.01, color = "#4C566A", key_glyph = draw_key_path)+
        geom_line(size = 1.2, key_glyph = draw_key_path) +
        scale_color_manual(values = cols) +
        #geom_smooth(method = "lm", se = FALSE) +
        facet_wrap(~hap, labeller=labeller(hap = hap_labels),
                   scales = "free_y") + 
        theme_simple(grid_lines = TRUE, axis_lines = TRUE) +
        scale_y_continuous(limits = c(0.05, 0.45)) +
        # scale_color_manual(values = c( "#5E81AC")) + #5E81AC "#94350b"
        ylab("Haplotype\nfrequency") +
        xlab("Birth cohort") +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.spacing = unit(2.5, "lines"),
              strip.text.x = element_blank(),
              panel.spacing.x = unit(0.5, "lines"),
              plot.margin = margin(b = 20, unit = "pt"))
# panel.grid.major.y = element_line(size = 0.1, color = "#4C566A")) 
p_freq


 # slopes
slopes <- map(gd_summaries, gd_slopes, n_founder_cohorts = 3,
              remove_founders = T) 

emp_slopes <- map_dfr(slopes, 2, .id = "hap") %>% 
                mutate(hap = factor(hap),
                hap = fct_relevel(hap, "chr5_6293", "chr7_12196","chr18_267"))

sim_slopes <- slopes %>% 
                map_dfr(1, .id = "hap") %>% 
                select(hap, Slope) %>% 
                mutate(hap = factor(hap),
                hap = fct_relevel(hap, "chr5_6293", "chr7_12196", "chr18_267"))

p_slopes <- sim_slopes %>% 
        ggplot(aes(Slope)) +
                geom_histogram(bins = 20, color = "black", size = 0.2,
                               fill = "#D8DEE9", alpha = 0.5) +
                theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
                scale_x_continuous(breaks = c(-0.005, 0, 0.005), 
                           labels = c("-0.005", "0", "0.005")) +
                scale_y_continuous(expand = c(0, 0)) +
                scale_fill_manual(values = cols) +
                facet_wrap(~hap, labeller=labeller(hap = hap_labels)) +
                geom_vline(data = emp_slopes, aes(xintercept = Slope, color = hap),
                           #color = "black",
                           size = 1.2, key_glyph = draw_key_path) +
                scale_color_manual(values = cols) +
                xlab("Slopes (directional selection)") +
                ylab("# simulations") +
                theme(strip.text.x = element_blank(),
                      #legend.position = "None",
                      plot.margin = margin(b = 20, unit = "pt")) 
              
                #ggtitle("Directional selection") +
                #theme(plot.title = element_text(size = 13, hjust = 0.5))
p_slopes   


# cumulative change
cumchange <- map(gd_summaries, gd_change, n_founder_cohorts = 3,
                 remove_founders = T) 
emp_change <- map_dfr(cumchange, 2, .id = "hap") %>% 
                mutate(hap = factor(hap),
                  hap = fct_relevel(hap, "chr5_6293", "chr7_12196","chr18_267"))

sim_change <- cumchange %>% 
                map_dfr(1, .id = "hap") %>% 
                select(hap, value) %>% 
                mutate(hap = factor(hap),
                  hap = fct_relevel(hap, "chr5_6293", "chr7_12196","chr18_267"))

p_change <- ggplot(sim_change, aes(value, fill = hap)) +
        geom_histogram(bins = 20, color = "black", size = 0.2,
                       fill = "#D8DEE9", alpha = 0.5) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1, 1.2)) + 
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = cols) +
        facet_wrap(~hap, labeller=labeller(hap = hap_labels)) +
        geom_vline(data = emp_change, aes(xintercept = value,color = hap),
                   #color = "black",
                   size = 1.2, key_glyph = draw_key_path) +
        scale_color_manual(values = cols) +
        xlab("Cumulative frequency change (balancing selection)") +
        ylab("# simulations") +
        theme(strip.text.x = element_blank(),
              #legend.position = "None",
              plot.margin = margin(b = 20, unit = "pt"))

library(patchwork)
p_freq / p_slopes / p_change +
        plot_layout(heights = c(2.5, 1, 1),
                    guides = 'collect') & theme(legend.position = 'bottom') &
        scale_color_viridis_d(labels=hap_labels, title="haplotype") 
        



plot_genedrop_cumulative_change(gd_summaries[[2]],
                        n_founder_cohorts = 3,
                        remove_founders = T)




get_slope <- function(sim) {
        sim <- sim %>% filter(birth_year > 1995)
        lm(p ~ birth_year, data = sim)$coefficients[[2]]
}

sim_slopes <- gds_df %>% 
        group_by(region, Simulation) %>% 
        nest() %>% 
        mutate(slope = map_dbl(data, get_slope)) %>% 
        ungroup()

emp_slopes <- gd_observed %>% 
                bind_rows(.id = "names") %>% 
                rename(birth_year = Cohort, region = names) %>% 
                group_by(region, Simulation) %>% 
                nest() %>% 
                mutate(slope = map_dbl(data, get_slope)) %>% 
                ungroup()
        
ggplot(sim_slopes, aes(slope)) +
        geom_histogram(bins = 50) +
        facet_wrap(~region) +
        geom_vline(data = emp_slopes, aes(xintercept = slope))

# cumulative change
out <- map(gd_summaries, gd_change) %>% 
        map_dfr(.[[2]], .id = "names") %>% 
        select(names, value)
out

dat <- gd_simulated[[1]] %>% 
        group_by(Simulation) %>% 
        group_map(~get_slope(.x)) %>% 
        unlist()

hist(dat)

out <- plot_genedrop_lm_slopes(gd_summary[[1]],
                        n_founder_cohorts = 5,
                        remove_founders = F)
out














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


