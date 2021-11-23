library(tidyverse)
# U = 2*genome_size*mutation_rate
make_slim(genome_size = 1e8, pop_size1 = 1000, pop_size2 = 200, 
          time1 = 10000, time2 = 11000,
          mut_rate_del = 1.2e-8, recomb_rate = 1e-8,
          mut1_dom_coeff = 0.2, mut1_gam_mean = -0.03, 
          mut1_gam_shape = 0.5, mut1_rel_freq = 0.95,
          mut2_dom_coeff = 0, mut2_val = -1, 
          mut2_rel_freq = 0.05,
          mut3_dom_coeff = NULL, mut3_gam_mean = NULL, 
          mut3_gam_shape = NULL, mut3_rel_freq = NULL,
          out_dir = "slim_sim/sims",
          seed = 100)

out_path = "slim_sim/sims"
seed = 100
# run slim
#system(paste0("slim -time -Memhist -seed ", seed, " slim_sim/sims/slim_code/sheep_", seed, ".slim"))
system(paste0("slim -time -seed ", seed," ", out_path, "/slim_code/sheep_", seed, ".slim"))

seed(101)
seeds <- sample(10, 100)


