library(tidyverse)
source("theme_simple.R")
muts <- read_delim("slim_sim/sims/muts/mutperind_100.txt")

test <- muts %>% 
        group_by(mut_id, pos, s, originG) %>% 
        tally() %>% 
        mutate(freq = n/400) %>% 
        ungroup() %>% 
        mutate(s_class = cut(s, breaks = c(0, -0.001, -0.01, c(seq(-0.1, -1.2, by = -0.1))))) %>% 
        group_by(s_class) %>% 
        summarise(mean(freq), min(freq), max(freq))
test

mut_classes <- muts %>% 
        group_by(mut_id, pos, s, originG) %>% 
        tally() %>% 
        mutate(freq = n/400) %>% 
        ungroup() %>% 
        mutate(s = ifelse(s < -1, -1, s)) %>% 
        mutate(s_class = cut(s, breaks = c(0, -0.001, -0.01, c(seq(-0.1, -1.1, by = -0.1)))),
               s_class = as.character(s_class)) %>% 
        rowwise() %>% 
        mutate(max_class = str_split_fixed(s_class, ",", 2)[1],
               max_class = as.numeric(str_remove(max_class, "\\("))) %>% 
        ungroup() %>% 
        mutate(max_class = ifelse(max_class == -1.1, -1, max_class))
        

ggplot(mut_classes, aes(s_class, freq)) +
        geom_jitter(width = 0.2, alpha = 0.3) 


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

pars <- expand_grid(freq = seq(0.01, 1, 0.01), s = c(0, -0.001, -0.01, c(seq(-0.1, -1, by = -0.1))))

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
        
        

