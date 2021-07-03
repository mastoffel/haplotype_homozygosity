# phase haploptypes
#library(tidyverse)
library(data.table)
library(glue)
library(furrr)
source("create_spec.R")

data_path <- "data"
out_path <- "output"
ped <- glue("{data_path}/ped.txt")

args_in <- commandArgs(trailingOnly=TRUE)
if (length(args_in) == 1) out_path <- args_in[1]

# peeling per chromosome
run_alpha_peel <- function(chr_num, data_path, out_path, ped) {
        
        # number of snps (-1 as data includes ID)
        num_snps <- ncol(fread(glue("{data_path}/genos/chr{chr_num}.txt"), nrows = 1))-1
        
        # make spec file
        genos <- glue("{data_path}/genos/chr{chr_num}.txt")
        spec_file <- glue("{data_path}/specs/spec_chr{chr_num}.txt")
        out_file <- glue("{out_path}/phased_hap_chr{chr_num}")
        
        create_spec(spec = spec_file,
                    nsnps = num_snps, 
                    ped = ped, 
                    genos = genos, 
                    out = out_file)
        
        # run alpha peel
        sys_command <- glue("./AlphaPeel_linux {spec_file}")
        system(sys_command)
}

plan(multiprocess, workers = 26)
future_walk(1:26, run_alpha_peel, data_path, out_path, ped)
