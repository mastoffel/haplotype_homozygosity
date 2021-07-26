library(tidyverse)
library(here)
source("theme_simple.R")
library(patchwork)
library(data.table)
library(GWASTools)
# read results from haplotype homozygosity scan
all_files <- list.files(here("output", "hap_results_imputed", "hap_len_200"), full.names = TRUE)
res_full <- map(all_files, read_delim, delim = "\t") %>% 
                bind_rows() %>% 
                rename(snp_num = snp_start)
qqPlot(res_full$p_val)
# any potential lethals?
res_full %>% 
        filter(obs == 0 & exp > 6)
res_full %>% 
        filter(chr == 18) %>% 
        arrange(p_val)

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
        filter(obs < exp) %>% 
        filter(p_val < 0.1) %>% 
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

cols <- viridis(2)
eff_tests <- 2*39149

p6 <- ggplot(gwas_plot, aes(positive_cum, -log10(p_val))) + 
        geom_hline(yintercept = -log10(0.05/(eff_tests)), linetype="dashed", color = "grey") +
        geom_point(data = gwas_plot %>% filter(-log10(p_val) <= -log10(0.05/(eff_tests))),
                   aes(color = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 0.8) +
        geom_point(data = gwas_plot %>% filter(-log10(p_val) > -log10(0.05/(eff_tests))), 
                   mapping = aes(fill = chromosome %%2 == 0), # aes(fill = direction),  0.00001
                   size = 2, shape = 21, stroke = 0.3, color = "grey") +
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,9), labels = as.character(0:8), breaks = 0:8) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) +
        scale_fill_manual(values = c("#ECEFF4","#d8dee9")) +
        scale_color_manual(values = c("#ECEFF4","#d8dee9")) + # #dbe1eb #d1d8e5  "#ECEFF4" #d8dee9
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        theme(axis.text = element_text(color = "black"), # axis.text.x size 8
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE, color = FALSE) +
        ggtitle("400K | Haplotype length: 400 SNPs")

p1


# run script a few times ...
p_final <- p1 / p2 / p3 / p4 / p5 / p6
ggsave("figs/manhattans_imputed.jpg", p_final, width = 6, height = 12)
