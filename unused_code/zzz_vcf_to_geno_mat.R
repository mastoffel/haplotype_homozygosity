library(vcfR)
library(purrr)
library(tidyverse)
vcf <- read.vcfR("data/sheep_phased_20.vcf", verbose = FALSE )
gt <- extract.gt(vcf)

gt

to_mat <- function(ind, gt) {
        genos <- strsplit(gt[, ind], "|", fixed = TRUE)
        genos2 <- unlist(genos)
        row1 <- as.numeric(genos2[c(TRUE, FALSE)])
        row2 <- as.numeric(genos2[c(FALSE, TRUE)])
        out <- rbind(row1, row2)
}

all_geno <- map(1:1000, to_mat, gt)

all_geno_mat <- matrix(unlist(all_geno), ncol = nrow(gt), byrow = TRUE)
dim(all_geno_mat)

# haplotypes

rowPaste <- function(row, cols, gt) {
        paste(gt[row, cols], collapse ="")
}

gt <- all_geno_mat
out <- map(1:nrow(gt), rowPaste, 6000:6020, gt)

sort(table(unlist(out)[1:2000]), decreasing = TRUE)

