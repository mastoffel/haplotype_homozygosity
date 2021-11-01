library(tidyverse)
library(here)
source("theme_simple.R")
library(patchwork)
library(data.table)
library(GWASTools)
library(ggeffects)
# read results from haplotype homozygosity scan
all_files <- list.files(here("output", "hap_results_imputed", "hap_len_500"), full.names = TRUE)

read_chrs <- function(file_path) {
        res <- read_delim(file_path, "\t") %>% 
                rename(snp_num = "snp_start") %>% 
                as_tibble() %>% 
                filter(obs < exp) %>% 
                # remove hap and give it simple id
                mutate(hap = 1:nrow(.)) %>% 
                filter(p_val < 0.05) 
        res
}

res_full <- map(all_files, read_chrs) %>% 
                bind_rows() 

#qqPlot(res_full$p_val)
# any potential lethals?
res_full %>% 
        filter(obs == 0 & exp > 9) %>% 
        print(n = 30)

# snp map
snp_map <- fread(here("data", "plink", "sheep.bim")) %>% 
        setNames(c("chr", "snp", "cM", "pos", "a", "b")) %>% 
        select(chr, snp, pos) %>% 
        group_by(chr) %>% 
        mutate(number = 1,
               snp_num = cumsum(number)) %>% 
        select(-number)
#
res <- res_full %>% 
        left_join(snp_map, by = c("chr", "snp_num"))

# chromosome info from assembly
chr_info <- read_delim("data/chromosome_info_oar31.txt", "\t") %>% 
        .[-1, ] %>% 
        rename(chromosome = Part) %>% 
        mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
        mutate(chromosome = as.integer(chromosome)) %>% 
        filter(!is.na(chromosome)) %>% 
        mutate(tot=cumsum(Length)-Length + cumsum(rep(35e6, 26))) %>% 
        dplyr::select(chromosome, tot) 

gwas_plot <- res %>% 
                rename(chromosome = chr) %>% 
                left_join(chr_info) %>% 
                arrange(chromosome, pos) %>%
                mutate(positive_cum = pos + tot) 

axisdf <- gwas_plot %>% 
                group_by(chromosome) %>% 
                summarise(center = (max(positive_cum) + min(positive_cum)) / 2 )

chr_labels <- c(c(1:10),"","12","", "14","", "16","", "18", "", "20","", "22","", "24","", "26")
chr_labels_full <- as.character(1:26)
cols <- c("#336B87", "#2A3132")

#cols <- viridis(2)
eff_tests <- 39149
gwas_plot <- gwas_plot %>% 
                group_by(chromosome) %>% 
                mutate(top_snp = ifelse(p_val == min(p_val), 1, 0)) %>% 
                group_by(chromosome, top_snp) %>% 
                mutate(top_snp2 = ifelse(top_snp == 1 & pos == min(pos) , 1, 0)) %>% 
                mutate(highlight = case_when(
                        top_snp2 == 1 & chromosome == 7 ~ 1,
                        top_snp2 == 1 & chromosome != 7 ~ 2,
                        top_snp2 != 1 & chromosome %% 2 == 0 ~ 3,
                        top_snp2 != 1 & chromosome %% 2 != 0 ~ 4
                )) %>% 
                mutate(highlight = factor(highlight))

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
        geom_point(data = gwas_plot %>% filter(highlight %in% c(1)) %>%
                           filter(p_val < 0.05/(eff_tests)),
                   size = 2.5, shape = 21, stroke = 0.1, fill = "#94350b", color = "black") +
        geom_point(data = gwas_plot %>% filter(highlight %in% c(2)) %>%
                           filter(p_val < 0.05/(eff_tests)),
                   size = 2.5, shape = 21, stroke = 0.1, fill = "#ccbe9b", color = "black") + 
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        theme(axis.text = element_text(color = "black"), # axis.text.x size 8
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE, color = FALSE)# +
# ggtitle("Haplotype length: 500 SNPs")

p_gwas

# p_gwas <- ggplot(gwas_plot, aes(positive_cum, -log10(p_val))) + 
#         geom_hline(yintercept = -log10(0.05/(eff_tests)), linetype="dashed", color = "grey") +
#         geom_point(data = gwas_plot %>% filter(-log10(p_val) <= -log10(0.05/(eff_tests))),
#                    aes(color = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
#                    size = 0.8) +
#         geom_point(data = gwas_plot %>% filter(-log10(p_val) > -log10(0.05/(eff_tests))), 
#                    mapping = aes(fill = chromosome %%2 == 0), # aes(fill = direction),  0.00001
#                    size = 2, shape = 21, stroke = 0.3, color = "grey") +
#         scale_x_continuous(labels = chr_labels, breaks= axisdf$center) +
#         scale_y_continuous(expand = c(0, 0), limits = c(0,9), labels = as.character(0:8), breaks = 0:8) +
#         xlab("Chromosome") + 
#         ylab(expression(-log[10](italic(p)))) +
#         scale_fill_manual(values = c("#ECEFF4","#d8dee9")) +
#         scale_color_manual(values = c("#ECEFF4","#d8dee9")) + # #dbe1eb #d1d8e5  "#ECEFF4" #d8dee9
#         theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
#         theme(axis.text = element_text(color = "black"), # axis.text.x size 8
#               axis.ticks = element_line(size = 0.1)) +
#         guides(fill=FALSE, color = FALSE)# +
       # ggtitle("Haplotype length: 500 SNPs")

#p_gwas


# run script a few times ...
#p_final <- p1 / p2 / p3 / p4 / p5 / p6 / p7
#ggsave("figs/manhattans_imputed2.jpg", p_final, width = 9, height = 14)


# haplotype frequency plots
haps_all <- read_delim("output/haps_all_500.txt")
p_freq <- haps_all %>% 
        mutate(region = factor(region),
               region = fct_reorder(region, chr)) %>% 
        filter(as.numeric(as.character(birth_year)) > 1990) %>% 
        group_by(birth_year, region) %>% 
        summarise(freq = sum(gt) / (n()*2)) %>% 
        ungroup() %>% 
        mutate(highlight = factor(ifelse(region == "chr7_12119", 1, 0))) %>% 
        ggplot(aes(birth_year, freq, group = 1, color = highlight)) +
        geom_line(size = 1.2) +
        #geom_smooth(method = "lm", se = FALSE) +
        facet_wrap(~region, ncol = 4) + 
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
        scale_color_manual(values = c("#ccbe9b", "#94350b")) +
        ylab("Haplotype\nfrequency") +
        xlab("Birth year") +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "none") 

# early failed pregnancies
failed_preg <- readRDS("output/ins_effects_500.rds")
marg_effs <- map(failed_preg, ggpredict, terms = c("mating_type")) 
plot_model(failed_preg[[4]], type = "pred", terms = c("mating_type"))
marg_df <- marg_effs %>% 
                map(as_tibble) %>% 
                bind_rows(.id = "region")
p_failed <- marg_df %>% 
        mutate(region = factor(region),
               region = fct_relevel(region, "chr18_267", after = 3)) %>% 
        mutate(mating_type = fct_recode(x, nn = "nonc_nonc", `cn|nc` = "c_nonc",
                                        cc = "c_c")) %>% 
        mutate(highlight = factor(ifelse(region == "chr7_12119", 1, 0))) %>% 
        #group_by(region) %>% 
        ggplot(aes(mating_type, predicted, color = highlight)) +
      
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 0.5)+
        geom_point(size = 2.5) +
        #geom_smooth(method = "lm", se = FALSE) +
        facet_wrap(~region, ncol = 4) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
        ylab("Probability of\nfailed insemination") +
        xlab("Mating type") +
       #scale_alpha_manual(values = c(0.5, 1)) +
        scale_color_manual(values = c("#ccbe9b", "#94350b")) +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "none") 
        # scale_color_manual(name="Mating type",
        #                    labels=c("nn (2 non-carriers)",
        #                             "cn|nc (1 carrier + 1 non-carrier)",
        #                             "cc (2 carriers"),
        #                    values=c("black","black","black"))
        

# survival plots
# surv_mods <- readRDS("output/surv_mods_hap500.rds")
# marg_effs <- ggpredict(surv_mods[[4]], terms = c("lamb", "gt"))
# plot_model(surv_mods[[1]], type = "pred", terms = c("lamb", "gt"))

p_final <- p_gwas / p_freq / p_failed   +
        plot_layout(heights = c(1.6, 1, 1)) +
        plot_annotation(tag_levels = "A")
ggsave("figs/haplotype_fig1.jpg", width = 7, height = 6)
