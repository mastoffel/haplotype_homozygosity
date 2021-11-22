library(tidyverse)

muts <- read_delim("output/sim_new.txt")

test <- muts %>% 
        group_by(mut_id, pos, s, originG) %>% 
        tally() %>% 
        mutate(freq = n/400) %>% 
        ungroup() %>% 
        mutate(s_class = cut(s, breaks = c(0, -0.001, -0.01, -0.1, -0.3, -0.5, -0.7, -1.2))) %>% 
        group_by(s_class) %>% 
        summarise(mean(freq))


make_slim(genome_size = 1e8, pop_size1 = 1000, pop_size2 = 200, 
          time1 = 10000, time2 = 11000,
          mut_rate_del = 7e-9, recomb_rate = 1e-8,
          mut1_dom_coeff = 0.2, mut1_gam_mean = -0.03, 
          mut1_gam_shape = 0.2,
          mut2_dom_coeff = 0, mut2_mean = -1, 
          mut2_sd = 0.05, mut2_rel_freq = 0.05,
          mut3_dom_coeff = NULL, mut3_gam_mean = NULL, 
          mut3_gam_shape = NULL, mut3_rel_freq = NULL,
          out_dir = "slim_sim/sims",
          seed = 100)
