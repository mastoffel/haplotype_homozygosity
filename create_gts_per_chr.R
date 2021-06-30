# make one genotype txt file per chr for alphapeel
#library(tidyverse)
library(glue)
library(data.table)
library(purrr)
# command line arguments
args_in <- commandArgs(trailingOnly=TRUE)
# plink files here
data_path <- "data"

if (length(args_in)==1) {
	out_path <- args_in[1]
} else if (length(args_in) == 0){
	out_path <- data_path
} else {
	stop("only one argument allowed (path to output")
}
	
# subfolder genos for single chromosome data
dir.create(file.path(out_path, "genos"))

# make alphapeel-ready txt files for genotypes per chromosome
create_gts_per_chr <- function(chr_num) {
        # extract genotype matrix
        system(glue("plink --bfile {data_path}/plink/sheep --recode A ",
                    "--chr {chr_num} --sheep --out {data_path}/genos/chr{chr_num}"))
        
        # read genotype file and modify for alpha peel
        genos <- fread(glue("{data_path}/genos/chr{chr_num}.raw"))
        # keep only ID and genotypes
        genos <- genos[, !c("FID", "PAT", "MAT", "SEX", "PHENOTYPE")]
        # change NA to 9 
        for (j in seq_len(ncol(genos))){
                set(genos,which(is.na(genos[[j]])),j,9)
        }
        # save as txt
        fwrite(genos, file = glue("{out_path}/genos/chr{chr_num}.txt"), 
               col.names = FALSE, sep = " ")  
        
}
# run
walk(1:26, create_gts_per_chr)
