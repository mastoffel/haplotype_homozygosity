library(tidyverse)
source("theme_simple.R")
library(ggdist)
library(ggrepel)

# 100 simulation runs 
files <- list.files("slim_sim/sims/muts/", full.names = TRUE)
full <- map_dfr(files, read_delim, .id = "run")

mut_classes <- full %>% 
        #group_by(run) %>% 
        #sample_frac(0.1) %>% 
        group_by(run, mut_id, pos, s, originG) %>% 
        tally() %>% 
        mutate(freq = n/400) %>% 
        ungroup() %>% 
        mutate(s_class = cut(s, breaks = c(0, -0.01, -0.05, -0.1, -0.2, -0.3, -0.4, -1.1))) %>% 
        #filter(!(s_class %in% c("(-0.8,-0.7]", "(-0.7,-0.6]", "(-0.6,-0.5]"))) %>% 
        mutate(s_class = fct_rev(s_class)) %>% 
        group_by(run, s_class) %>% 
        summarise(mf = mean(freq), lf = quantile(freq, 0.01), 
                  hf = quantile(freq, 0.9), 
                  hf2 = quantile(freq, 0.99),
                                 n = n()) %>% 
        group_by(s_class) %>% 
        summarise(mf = mean(mf), hf = mean(hf), hf2 = mean(hf2), lf = mean(lf), n = mean(n)) %>%  
        mutate(n = round(n, 0)) 
        

pd <- position_dodge(0.1)

mut_classes %>% 
        ggplot(aes(x = mf, y = s_class, size = n)) + 
        geom_point() +
       # geom_point(aes(x = hf, y = s_class), 
        #           size = 2, shape = "|", stroke = 1) +
        geom_linerange(aes(xmin = mf, xmax = hf), size = 1.5) +
        geom_linerange(aes(xmin = mf, xmax = hf2), size = 0.5) +
        theme_simple(grid_lines = FALSE) +
        geom_label(label = mut_classes$n,
                   nudge_x = 0.1, nudge_y = 0.4, 
                    size = 3) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        xlab("Mutation frequency (Mean, 90th, 99th percentile)") +
        theme(legend.position = "none")
        #xlim(c(0, 0.5))
        
        


p <- ggplot(mut_classes, aes(x = freq, y=s_class)) +
        stat_interval(
                     position = position_nudge(y = -0.2),
                     orientation = "horizontal",
                     .width = c(.5,.8,.99),
                     stroke = 2,
                     size = c(6)) +
        theme_simple(grid_lines = FALSE) +
        xlim(c(0, 0.3)) +
        scale_colour_viridis_d("Percentile", direction = -1) +
        stat_dots(size = 2,
                  scale = 0.65,
                  color = "#2E3440",
                  fill = "#4C566A",
                  shape=23,
                  stroke = 0.001,
                  position = "dodge",
                  alpha = 0.1)
p
        
ggsave("figs/sims.jpg",p, width = 6, height = 6)

stat_dots(position = "dodge", orientation = "horizontal", 
                  layout = "swarm",
                  scale = 0.6,
                  size = 3,
                  color = "#6f4e37",
                  fill = "#6f4e37")
        
#stat_slab()
        stat_dots(
                quantiles = NA,
                orientation = "horizontal",
                #normalize = "none",
                #scale = .87,
                scale = 1,
                dotsize = 0.01,
                stackratio = 0.1,
               # layout = "swarm",
                color = "#6f4e37",
                fill = "#6f4e37"
        )




mut_classes %>% 
        group_by(run, s_class) %>% 
        summarise(mean_freq = mean(freq), min_freq = min(freq), max_freq=max(freq),
                  mean_G = mean(originG), mean_s = mean(s)) %>% 
        rowwise() %>% 
        mutate(max_class = str_split_fixed(s_class, ",", 2)[1],
               max_class = as.numeric(str_remove(max_class, "\\("))) %>% 
        mutate(max_class = ifelse(max_class == -1.1, -1, max_class))

mut_classes2 <- mut_classes %>% 
        group_by(s_class) %>% 
        summarise(mean_freq = mean(mean_freq), min_freq = mean(min_freq),max_freq = mean(max_freq))

ggplot(mut_classes2, aes(s_class, mean_freq)) +
        geom_pointrange(aes(ymin = min_freq, ymax = max_freq), size = 0.3) + 
        scale_y_log10()


# expected
s = 0.3
freq = 0.2
E_k_hom <- 500
O_k_hom <- (1-s) * 500
n <- 5952
O_k_nonhom <- n-O_k_hom
E_k_nonhom <- n-E_k_hom

# 5952
chisq <- ((O_k_hom - E_k_hom)^2)/E_k_hom + ((O_k_nonhom - E_k_nonhom)^2)/E_k_nonhom
chi_p <- pchisq(chisq, df = 1, lower.tail = FALSE)
chi_p 

get_p <- function(freq, s, n = 5952) {
        E_k_hom <- freq^2 * n
        O_k_hom <- (1-s) * E_k_hom
        O_k_nonhom <- n-O_k_hom
        E_k_nonhom <- n-E_k_hom
        chisq <- ((O_k_hom - E_k_hom)^2)/E_k_hom + ((O_k_nonhom - E_k_nonhom)^2)/E_k_nonhom
        chi_p <- pchisq(chisq, df = 1, lower.tail = FALSE)
        out <- chi_p * 39149
        out <- ifelse(out>1, 1, out)
}

get_p(freq = 0.25, -0.5)

pars <- expand_grid(freq = seq(0.01, 1, 0.01), s = c(0, -0.0001, -0.001, -0.01, c(seq(-0.1, -1, by = -0.1))))

pars$p <- map2(pars$freq,pars$s, get_p) %>% 
        unlist()

pars_thres <- pars %>% 
        group_by(s) %>% 
        filter(p <= 0.05) %>% 
        filter(freq == min(freq)) %>% 
        ungroup()
pars_thres[8, 2] <- -0.3

df <- mut_classes %>% 
        rename(mut_freq = freq) %>% 
        left_join(pars_thres, by = c("max_class" = "s"))
df_sum <- df %>% 
                group_by(s_class) %>% 
                summarise(threshold = mean(freq))

ggplot(df, aes(s_class, mut_freq)) +
        geom_jitter(width = 0.2, alpha = 0.3) +
        geom_point(aes(s_class, threshold), data = df_sum, color = "blue", size = 2) +
        theme_simple()
        
        

