# phase haploptypes
library(tidyverse)
library(data.table)
library(glue)
source("create_spec.R")

data_path <- "data"
out_path <- "output"
ped <- glue("{data_path}/ped.txt")



# peeling per chromosome
chr_num <- 24

# number of snps (-1 as data includes ID)
num_snps <- ncol(fread(glue("{data_path}/genos/chr{chr_num}.txt"), nrows = 1))-1

# make spec file
genos <- glue("{data_path}/genos/chr{chr_num}.txt")
create_spec(spec = glue("{data_path}/specs/spec_chr{chr_num}.txt"),
            nsnps = num_snps, 
            ped = ped, 
            genos = genos, 
            out = glue("{out_path}/phased_hap_chr{chr_num}"))

