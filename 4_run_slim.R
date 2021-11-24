library(tidyverse)
library(glue)
source("make_slim.R")
library(furrr)

args_in <- commandArgs(trailingOnly=TRUE)
out_path <- args_in[1]

# U = 2*genome_size*mutation_rate
set.seed(111)
seeds <- sample(1:10000, 100, replace = FALSE)

run_slim <- function(seed) {
        
        make_slim(genome_size = 1e2, pop_size1 = 1000, pop_size2 = 200, 
                  time1 = 10000, time2 = 11000,
                  mut_rate_del = 1.2e-8, recomb_rate = 1e-8,
                  mut1_dom_coeff = 0.2, mut1_gam_mean = -0.03, 
                  mut1_gam_shape = 0.5, mut1_rel_freq = 0.95,
                  mut2_dom_coeff = 0, mut2_val = -1, 
                  mut2_rel_freq = 0.05,
                  mut3_dom_coeff = NULL, mut3_gam_mean = NULL, 
                  mut3_gam_shape = NULL, mut3_rel_freq = NULL,
                  out_dir = glue("{out_path}/slim_sim/sims"),
                  seed = seed)
        
        system(glue("slim -time -seed {seed} {out_path}/slim_sim/sims/slim_code/sheep_{seed}.slim"))
}

plan(multisession, workers = 2)
walk(seeds[1:4], run_slim)
