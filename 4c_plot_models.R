library(tidyverse)
library(brms)
library(broom.mixed)
library(here)
library(ggdist)
library(patchwork)
source("theme_simple.R")

# survival plots
surv <- read_delim(here("output", "survival_marginal_effects.txt"))

# get post. mean + CI
surv %>%
  group_by(predictor) %>%
  summarise(post_mean = mean(estimate),
            CI_lower = quantile(estimate, probs = 0.025),
            CI_upper = quantile(estimate, probs = 0.975))

cols <- c('#161032', '#C73E1D', "#FFC53A")
hap_labels <- c("SEL05", "SEL07", "SEL18")

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
  scale_color_manual("Haplotype", labels = hap_labels, values = cols) +
  scale_fill_manual("Haplotype", labels = hap_labels,values = cols) +
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
        #legend.position = "none",
        panel.spacing = unit(1, "lines"),
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
  scale_color_manual("Haplotype", labels = hap_labels, values = cols) +
  scale_fill_manual("Haplotype", labels = hap_labels,values = cols) +
  facet_grid(~chr ,scales = "free_x") +
  scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-0.5, 0.8)) +
  ylab("Haplotype\ncopies") +
  xlab("Predicted change in body weight (kg)") +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        #legend.position = "none",
        panel.spacing = unit(1, "lines"),
        strip.text.x = element_blank())
p_weight

# final plot
p_final <-  p_surv / p_weight  +
  plot_annotation(tag_levels = "A")  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
        #text = element_text("sans"))
p_final

ggsave("figs/models.tiff", width = 5.4, height = 4, dpi = 600)
ggsave("figs/models_pub.pdf", width = 5.4, height = 4, device = cairo_pdf)



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
