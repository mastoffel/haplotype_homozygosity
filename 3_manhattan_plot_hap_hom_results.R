library(tidyverse)
library(here)
source("theme_simple.R")
library(patchwork)
library(data.table)
library(GWASTools)
library(ggeffects)
library(tidybayes)
library(genedroppeR)
# read results from haplotype homozygosity scan
all_files <- list.files(here("output", "hap_results_imputed", "hap_len_400"), full.names = TRUE)

read_chrs <- function(file_path) {
        res <- read_delim(file_path, "\t") %>% 
                rename(snp_num = "snp_start") %>% 
                as_tibble() %>% 
                filter(obs < exp) %>% 
                # remove hap and give it simple id
                mutate(hap = 1:nrow(.)) #%>% 
                #filter(p_val < 0.05) 
        res
}

res_full <- map(all_files, read_chrs) %>% 
                bind_rows() 

#qqPlot(res_full$p_val)

# snp map
snp_map <- fread(here("data", "plink", "sheep.bim")) %>% 
        setNames(c("chr", "snp", "cM", "pos", "a", "b")) %>% 
        select(chr, snp, pos) %>% 
        group_by(chr) %>% 
        mutate(number = 1,
               snp_num = cumsum(number)) %>% 
        select(-number)

# any potential lethals?
res_full %>% 
  filter(obs == 0 & exp > 8) %>% 
  print(n = 30)

# start and end
snp_map %>% filter(chr == 7, snp_num == 1072 | snp_num == 1072+399) %>% mutate(pos_mb = pos/1e6)

#
res <- res_full %>% 
        left_join(snp_map, by = c("chr", "snp_num"))

# plots
cols <- c("#33658A", "#86BBD8", "#2F4858")
#c( "#D08770",  "#5E81AC","#A3BE8C")

# chromosome info from assembly
chr_info <- read_delim("data/chromosome_info_oar31.txt", "\t") %>% 
        .[-1, ] %>% 
        rename(chromosome = Part) %>% 
        mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
        mutate(chromosome = as.integer(chromosome)) %>% 
        filter(!is.na(chromosome)) %>% 
        mutate(tot=cumsum(Length)-Length + cumsum(rep(35e6, 26))) %>% 
        dplyr::select(chromosome, tot) 

gwas_plot_tmp <- res %>% 
                rename(chromosome = chr) %>% 
                left_join(chr_info) %>% 
                arrange(chromosome, pos) %>%
                mutate(positive_cum = pos + tot) 

axisdf <- gwas_plot_tmp %>% 
                group_by(chromosome) %>% 
                summarise(center = (max(positive_cum) + min(positive_cum)) / 2 )

chr_labels <- c(c(1:10),"","12","", "14","", "16","", "18", "", "20","", "22","", "24","", "26")
chr_labels_full <- as.character(1:26)
#cols <- c("#336B87", "#2A3132")

#cols <- viridis(2)
eff_tests <- 39149 #39149
gwas_plot <- gwas_plot_tmp %>% 
                #filter(p_val < 0.01) %>% 
                group_by(chromosome) %>% 
                mutate(top_snp = ifelse(p_val == min(p_val), 1, 0)) %>% 
                group_by(chromosome, top_snp) %>% 
                mutate(top_snp2 = ifelse(top_snp == 1 & pos == min(pos) , 1, 0)) 
  
library(ggnewscale)
p_gwas <- ggplot(gwas_plot, aes(positive_cum, -log10(p_val))) + 
        geom_hline(yintercept = -log10(0.05/(eff_tests)), linetype="dashed", color = "grey") +
        geom_point(data = gwas_plot %>% filter(-log10(p_val) <= -log10(0.05/(eff_tests))),
                   aes(color = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 0.8) +
        geom_point(data = gwas_plot %>% filter(-log10(p_val) > -log10(0.05/(eff_tests))), 
                   mapping = aes(fill = chromosome %%2 == 0), # aes(fill = direction),  0.00001
                   size = 2, shape = 21, stroke = 0.3, color = "grey") +
        # "#ccbe9b", "#94350b"
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,9), labels = as.character(0:8), breaks = 0:8) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) +
        scale_fill_manual(values = c("#ECEFF4","#d8dee9")) +
        scale_color_manual(values = c("#ECEFF4","#d8dee9")) + # #dbe1eb #d1d8e5  "#ECEFF4" #d8dee9
        guides(fill=FALSE, color = FALSE) +
        new_scale_fill() +
        geom_point(data = gwas_plot %>% 
                           filter(p_val < 0.05/(eff_tests)) %>% 
                           filter(top_snp2 == 1),
                   size = 3, shape = 21, stroke = 0.1, mapping = aes(fill = snp)) + # "#94350b"
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        scale_fill_manual(values = cols[c(3,1,2)]) +
        theme(axis.text = element_text(color = "black"), # axis.text.x size 8
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE, color = FALSE)# +
# ggtitle("Haplotype length: 500 SNPs")
#p_gwas


hap_names <- c(
  "chr5_6293" = "SEL05",
  "chr7_12196" = "SEL07",
  "chr18_267" = "SEL18"
)

# row with barplots?
top_haps <- read_delim(here("output", "top_haps_400.txt")) %>% 
              select(region, obs, exp) %>% 
              mutate(prop_obs = paste0("-", round((1-obs/exp)*100, 0), "%")) %>% 
              pivot_longer(cols = obs:exp,
                           names_to = "type",
                           values_to = "ind_count") %>% 
              mutate(prop_obs = ifelse(type == "exp", "", paste0("(", prop_obs, ")"))) %>% 
              mutate(count_label = paste0(round(ind_count,0))) %>% 
              mutate(region = factor(region, levels = c("chr5_6293", "chr7_12196",
                                                        "chr18_267"))) %>% 
              mutate(type_region = paste0(region, type))

p_num <- ggplot(top_haps, aes(type, ind_count, col = type_region)) +
  #geom_col(fill = "#ccbe9b") +
  geom_point(size = 2.5) +
  geom_segment(aes(x=type, xend=type, y=0, yend=ind_count)) +
  facet_wrap(~region, labeller=labeller(region = hap_names)) +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
  scale_y_continuous(limits = c(0, 270)) +
 # scale_color_manual(values = c("#ccbe9b", "#94350b")) +
  #ylab("Haplotype\nfrequency") +
  ylab("# Homozygous\noffspring") +
  xlab("Category") +
  scale_color_manual(values = c("#2E3440", cols[3], "#2E3440", cols[1],
                                "#2E3440", cols[2])) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2.5, "lines")) +
  geom_text(aes(label = count_label), vjust = 0.4, hjust = 1.5,
            size = 3) +
  geom_text(aes(label = prop_obs), vjust = 1.8, hjust = 1.1,
            size = 3)
p_num
# run script a few times ...
#p_final <- p1 / p2 / p3 / p4 / p5 / p6 / p7
#ggsave("figs/manhattans_imputed2.jpg", p_final, width = 9, height = 14)


# haplotype frequency plots
haps_all <- read_delim(here("output", "haps400_and_fitness.txt"))

# get genedrops
hap_names_full <- c("chr18_267", "chr7_12196",
               "chr5_6293")
load_gd <- function(hap_name) {
  gd <- read_delim(here("output", paste0("genedrop_hap_", hap_name, ".txt")))
  gd_summ <- summary_genedrop(gd)[[2]] %>% 
              as_tibble() %>% 
              mutate(region = hap_name)
}
genedrops <- map_dfr(hap_names_full, load_gd) %>% 
              rename(birth_year = Cohort, freq = Count) %>% 
              mutate(region = factor(region, levels = c("chr5_6293", "chr7_12196",
                                            "chr18_267")))

p_freq <- haps_all %>% 
        mutate(region = factor(region),
               region = fct_reorder(region, chr)) %>% 
        filter(as.numeric(as.character(birth_year)) >= 1990) %>% 
        group_by(birth_year, region) %>% 
        summarise(freq = sum(gt) / (n()*2)) %>% 
        ungroup() %>% 
        #mutate(highlight = factor(ifelse(region == "chr7_12119", 1, 0))) %>% 
        ggplot(aes(birth_year, freq, group = 1, color = region)) +
        geom_line(data=genedrops, alpha = 0.7, aes(y = p, group = Simulation),
                  size = 0.05, color = "#D8DEE9")+
        geom_line(size = 1.2) +
        scale_color_manual(values = cols) +
        #geom_smooth(method = "lm", se = FALSE) +
        facet_grid(~region, labeller=labeller(region = hap_names)) + 
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
       # scale_color_manual(values = c( "#5E81AC")) + #5E81AC "#94350b"
        ylab("Haplotype\nfrequency") +
        xlab("Birth year") +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "none",
              panel.spacing = unit(2.5, "lines"),
              strip.text.x = element_blank())
             # panel.grid.major.y = element_line(size = 0.1, color = "#4C566A")) 
p_freq
#ggsave("figs/hap_freq.jpg", p_freq, width = 7, height = 2)


# survival plots
surv <- read_delim(here("output", "survival_marginal_effects.txt"))

# get post. mean + CI
surv %>% 
  group_by(predictor) %>% 
  summarise(post_mean = mean(estimate),
            CI_lower = quantile(estimate, probs = 0.025),
            CI_upper = quantile(estimate, probs = 0.975))

p_surv <- surv %>% 
  mutate(chr = str_remove(hap, "gt"))%>% 
  mutate(chr = case_when(
    chr == "5" ~ "hap05",
    chr == "7" ~ "hap07",
    chr == "18" ~ "hap18"
  )) %>% 
  mutate(chr = factor(chr, levels = c("hap05", "hap07", "hap18"))) %>% 
  mutate(copies = str_sub(predictor, start = -1)) %>% 
  ggplot(aes(x = estimate, y = copies,
    color = chr,
    fill = chr)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 1) +
  stat_halfeye(#mapping = aes(fill=stat(
    #cut_cdf_qi(cdf, .width = c(.66,.95,1)))),
    #interval_color = "#4C566A",
    #point_color = "#4C566A",
    #slab_color = "#4C566A",
    adjust = 3,
    #width = .6,
    #.width = 0, 
    justification = -.1, 
    height = 0.8,
    slab_size = 0.5,
    slab_alpha = 0.7
    #slab_fill = "#E5E9F0"
  ) +

  #scale_slab_color_discrete(cols) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
 # scale_slab_color_discrete(values = cols) +
  #scale_fill_manual(values=cols)+
  #scale_slab_color_discrete() +
  #scale_fill_manual(values = c("#D08770", "#5E81AC","#A3BE8C")) +
  facet_grid(~chr ,scales = "free_x") +
  scale_x_continuous(breaks = seq(-30, 30, 10), limits = c(-31, 31), 
                     labels = c("", "-20%", "", "0%", "", "20%", "")) + #limits = c(-36, 36)
  ylab("Haplotype\ncopies") +
  xlab("Predicted change in first-year survival") +
  #scale_fill_brewer(direction = -1, na.translate = FALSE) +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2.5, "lines"),
        strip.text.x = element_blank()) 
p_surv


# weight model plots
fit <- readRDS("output/haps_weight_mod.RDS")
post_df <- posterior_samples(fit) %>% 
  as_tibble() %>% 
  select(b_gt181:b_gt72) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(#copies = str_sub(name, -1),
         hap = str_sub(name, start=3, end=-2)) %>% 
  mutate(chr = str_remove(hap, "gt"))%>% 
  mutate(chr = case_when(
    chr == "5" ~ "hap05",
    chr == "7" ~ "hap07",
    chr == "18" ~ "hap18"
  )) %>% 
  rename(predictor = name,
         estimate = value) %>% 
  mutate(chr = factor(chr, levels = c("hap05", "hap07", "hap18"))) %>% 
  mutate(copies = str_sub(predictor, start = -1))


p_weight <- post_df %>% 
  ggplot(aes(x = estimate, y = copies,
           color = chr,
           fill = chr)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 1) +
  stat_halfeye(
    adjust = 3,
    justification = -.1, 
    height = 0.8,
    slab_size = 0.5,
    slab_alpha = 0.7
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  facet_grid(~chr ,scales = "free_x") +
  scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-0.5, 0.8)) +
  ylab("Haplotype\ncopies") +
  xlab("Predicted change in body weight (kg)") +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.spacing = unit(2.5, "lines"),
        strip.text.x = element_blank()) 
p_weight
# final plot
p_final <- p_gwas / p_num/ p_freq / p_weight / p_surv   +
  plot_layout(heights = c(2, 1.3, 1.5, 1.5, 1.5)) +
  plot_annotation(tag_levels = "A")
#p_final

ggsave("figs/haplotype_fig2.jpg", width = 6.5, height = 9)




# survival plots
surv_mods <- readRDS("output/firstyear_surv_mods_hap500.rds")
map(surv_mods, tidy, conf.int = TRUE)
#marg_effs <- map(surv_mods, ggpredict, terms = c("gt")) 
marg_effs <- map(surv_mods, ggemmeans, terms = c("gt")) 
marg_df <- marg_effs %>% 
  map(as_tibble) %>% 
  bind_rows(.id = "region")

p_surv <- marg_df %>% 
  mutate(region = factor(region),
         region = fct_relevel(region, "chr18_267", after = 3)) %>% 
  #mutate(mating_type = fct_recode(x, nn = "nonc_nonc", `cn|nc` = "c_nonc",
  #                                cc = "c_c")) %>% 
  #mutate(highlight = factor(ifelse(region == "chr9_6571", 1, 0))) %>% 
  #group_by(region) %>% 
  ggplot(aes(x, predicted, color = highlight)) +
  
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 0.5)+
  geom_point(size = 2.5) +
  #geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~region, ncol = 4) +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
  ylab("Survival\nprobability") +
  xlab("Genotype (0=hom, 1=het, 2=hom)") +
  #scale_alpha_manual(values = c(0.5, 1)) +
  scale_color_manual(values = c("#ccbe9b", "#94350b")) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.spacing = unit(1.5, "lines")) 

p_final <- p_gwas / p_freq / p_surv   +
  plot_layout(heights = c(1.6, 1, 1)) +
  plot_annotation(tag_levels = "A")
p_final
#plot_model(surv_mods[[1]], type = "pred", terms = "gt")
ggsave("figs/haplotype_fig1_nolabs.jpg", width = 7, height = 6)
