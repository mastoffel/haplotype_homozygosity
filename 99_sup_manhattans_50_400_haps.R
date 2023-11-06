# make manhattan plots for haplotype length 50-400

library(tidyverse)
library(here)
source("theme_simple.R")
library(patchwork)
library(data.table)

# snp_map
snp_map <- fread(here("data", "plink", "sheep.bim")) %>% 
        setNames(c("chr", "snp", "cM", "pos", "a", "b")) %>% 
        select(chr, snp, pos) %>% 
        group_by(chr) %>% 
        mutate(number = 1,
               snp_num = cumsum(number)) %>% 
        select(-number)


get_plots <- function(hap_length, snp_map) {
        
        # read results from haplotype homozygosity scan
        all_files <- list.files(here("output", "hap_results_imputed", paste0("hap_len_", hap_length)), 
                                full.names = TRUE)
        
        read_chrs <- function(file_path) {
                res <- read_delim(file_path, "\t") %>% 
                        rename(snp_num = "snp_start") %>% 
                        as_tibble() %>% 
                        filter(obs < exp) %>% 
                        mutate(hap = 1:nrow(.)) #%>% 
                #filter(p_val < 0.01) 
                res
        }
        
        res_full <- map(all_files, read_chrs) %>% 
                bind_rows() 
        
        # add snp map
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
        eff_tests <- 39149 #39149
        
        gwas_plot <- gwas_plot_tmp %>% 
                #filter(p_val < 0.01) %>% 
                group_by(chromosome) %>% 
                mutate(top_snp = ifelse(p_val == min(p_val), 1, 0)) %>% 
                group_by(chromosome, top_snp) %>% 
                mutate(top_snp2 = ifelse(top_snp == 1 & pos == min(pos) , 1, 0)) 
        
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
                theme_simple(grid_lines = FALSE) +
                theme(axis.text = element_text(color = "black"), # axis.text.x size 8
                      axis.ticks = element_line(size = 0.1)) +
                guides(fill=FALSE, color = FALSE)# +
        # ggtitle("Haplotype length: 500 SNPs")
        p_gwas
        
}

hap_lengths <- c(100, 200, 300, 400, 500)
all_plots <- map(hap_lengths, get_plots, snp_map)
p <- all_plots
p_final <- p[[1]] / p[[2]] / p[[3]] / p[[4]] / p[[5]]
ggsave(here("figs/manhattans_hap_length.jpg"), p_final, width = 7.5, height = 12)
       